//! Merge multiple fastp JSON output files into a single aggregated result,
//! as if fastp had been run on the entire dataset at once.

use serde_json::{Map, Value};

/// Sum arrays element-wise, padding shorter arrays with zero.
///
/// Preserves integer values when all inputs are integers; uses `f64` otherwise.
#[must_use]
pub fn sum_arrays(arrays: &[&[Value]]) -> Vec<Value> {
    if arrays.is_empty() {
        return Vec::new();
    }
    let max_len = arrays.iter().map(|a| a.len()).max().unwrap_or(0);

    let all_integers = arrays
        .iter()
        .all(|arr| arr.iter().all(|v| v.is_i64() || v.is_u64()));

    if all_integers {
        let mut result = vec![0_i64; max_len];
        for arr in arrays {
            for (i, val) in arr.iter().enumerate() {
                result[i] = result[i].saturating_add(val.as_i64().unwrap_or(0));
            }
        }
        result.into_iter().map(Value::from).collect()
    } else {
        let mut result = vec![0_f64; max_len];
        for arr in arrays {
            for (i, val) in arr.iter().enumerate() {
                result[i] += val.as_f64().unwrap_or(0.0);
            }
        }
        result.into_iter().map(Value::from).collect()
    }
}

/// Average arrays element-wise, handling variable-length inputs gracefully.
#[must_use]
#[allow(clippy::cast_precision_loss)]
pub fn average_arrays(arrays: &[&[Value]]) -> Vec<Value> {
    if arrays.is_empty() {
        return Vec::new();
    }
    let max_len = arrays.iter().map(|a| a.len()).max().unwrap_or(0);
    let mut sums = vec![0_f64; max_len];
    let mut counts = vec![0_usize; max_len];
    for arr in arrays {
        for (i, val) in arr.iter().enumerate() {
            sums[i] += val.as_f64().unwrap_or(0.0);
            counts[i] += 1;
        }
    }
    sums.iter()
        .zip(counts.iter())
        .map(|(&sum, &count)| {
            if count > 0 {
                Value::from(sum / count as f64)
            } else {
                Value::from(0.0)
            }
        })
        .collect()
}

/// Merge kmer count objects by summing counts per kmer.
#[must_use]
pub fn merge_kmer_counts(dicts: &[&Map<String, Value>]) -> Map<String, Value> {
    let mut merged: Map<String, Value> = Map::new();
    for dict in dicts {
        for (kmer, count) in *dict {
            let existing = merged.get(kmer).and_then(Value::as_i64).unwrap_or(0);
            merged.insert(
                kmer.clone(),
                Value::from(existing + count.as_i64().unwrap_or(0)),
            );
        }
    }
    merged
}

/// Merge quality or content curve objects by averaging each key's array across chunks.
#[must_use]
pub fn merge_curve_maps(curves: &[&Map<String, Value>]) -> Map<String, Value> {
    if curves.is_empty() {
        return Map::new();
    }
    let mut merged = Map::new();
    for key in curves[0].keys() {
        let arrays: Vec<&[Value]> = curves
            .iter()
            .filter_map(|c| c.get(key)?.as_array().map(Vec::as_slice))
            .collect();
        merged.insert(key.clone(), Value::Array(average_arrays(&arrays)));
    }
    merged
}

/// Merge adapter count objects by summing counts per adapter sequence.
#[must_use]
pub fn merge_adapter_counts(dicts: &[&Map<String, Value>]) -> Map<String, Value> {
    let mut merged: Map<String, Value> = Map::new();
    for dict in dicts {
        for (adapter, count) in *dict {
            let existing = merged.get(adapter).and_then(Value::as_i64).unwrap_or(0);
            merged.insert(
                adapter.clone(),
                Value::from(existing + count.as_i64().unwrap_or(0)),
            );
        }
    }
    merged
}

/// Merge filtering result objects by summing all count fields.
#[must_use]
pub fn merge_filtering_result(results: &[Value]) -> Value {
    const FIELDS: &[&str] = &[
        "passed_filter_reads",
        "low_quality_reads",
        "too_many_N_reads",
        "too_short_reads",
        "too_long_reads",
    ];
    let mut merged = Map::new();
    for field in FIELDS {
        let sum: i64 = results
            .iter()
            .filter_map(|r| r.get(*field)?.as_i64())
            .sum();
        merged.insert((*field).to_string(), Value::from(sum));
    }
    Value::Object(merged)
}

/// Merge insert size objects by summing histograms and recomputing the peak.
#[must_use]
pub fn merge_insert_size(sizes: &[Value]) -> Value {
    let histograms: Vec<&[Value]> = sizes
        .iter()
        .filter_map(|s| s.get("histogram")?.as_array().map(Vec::as_slice))
        .collect();
    let merged_histogram = sum_arrays(&histograms);

    let unknown: i64 = sizes
        .iter()
        .filter_map(|s| s.get("unknown")?.as_i64())
        .sum();

    let peak = merged_histogram
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| {
            a.as_f64()
                .partial_cmp(&b.as_f64())
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map_or(0, |(i, _)| i);

    let mut merged = Map::new();
    merged.insert("peak".to_string(), Value::from(peak));
    merged.insert("unknown".to_string(), Value::from(unknown));
    merged.insert("histogram".to_string(), Value::Array(merged_histogram));
    Value::Object(merged)
}

/// Compute a weighted-average duplication rate, weighting each chunk by its total read count.
#[must_use]
pub fn merge_duplication(dups: &[Value], summaries: &[Value]) -> Value {
    let mut weighted_sum = 0.0_f64;
    let mut total_reads = 0.0_f64;
    for (dup, summary) in dups.iter().zip(summaries.iter()) {
        let reads = summary
            .get("before_filtering")
            .and_then(|b| b.get("total_reads"))
            .and_then(Value::as_f64)
            .unwrap_or(0.0);
        let rate = dup.get("rate").and_then(Value::as_f64).unwrap_or(0.0);
        weighted_sum += rate * reads;
        total_reads += reads;
    }
    let rate = if total_reads > 0.0 {
        weighted_sum / total_reads
    } else {
        0.0
    };
    let mut merged = Map::new();
    merged.insert("rate".to_string(), Value::from(rate));
    Value::Object(merged)
}

/// Recalculate derived rate and mean-length fields for a summary section after totals are merged.
/// Only called for `before_filtering` / `after_filtering` sections, not per-read sections.
#[allow(clippy::cast_possible_truncation)]
fn apply_summary_rates(merged: &mut Map<String, Value>, stats: &[Value]) {
    let total_bases = merged
        .get("total_bases")
        .and_then(Value::as_f64)
        .unwrap_or(0.0);
    let total_reads = merged
        .get("total_reads")
        .and_then(Value::as_f64)
        .unwrap_or(0.0);

    if total_bases > 0.0 {
        let q20 = merged.get("q20_bases").and_then(Value::as_f64).unwrap_or(0.0);
        let q30 = merged.get("q30_bases").and_then(Value::as_f64).unwrap_or(0.0);
        merged.insert("q20_rate".to_string(), Value::from(q20 / total_bases));
        merged.insert("q30_rate".to_string(), Value::from(q30 / total_bases));
    } else {
        merged.insert("q20_rate".to_string(), Value::from(0.0));
        merged.insert("q30_rate".to_string(), Value::from(0.0));
    }

    // Weighted mean read lengths. fastp reports `total_reads` counting both mates
    // for paired-end data, but that factor is identical across chunks and cancels
    // in the weighted average, so weighting directly by `total_reads` is correct
    // for both single- and paired-end. `read2_mean_length` is absent for single-end
    // summaries, so each field is emitted only when the input actually contains it.
    if total_reads > 0.0 {
        for field in &["read1_mean_length", "read2_mean_length"] {
            if stats.first().and_then(|s| s.get(*field)).is_none() {
                continue;
            }
            let weighted_sum: f64 = stats
                .iter()
                .filter_map(|s| {
                    let mean = s.get(*field)?.as_f64()?;
                    let reads = s.get("total_reads")?.as_f64()?;
                    Some(mean * reads)
                })
                .sum();
            // Truncation is intentional: fastp stores mean lengths as integers
            merged.insert(
                (*field).to_string(),
                Value::from((weighted_sum / total_reads) as i64),
            );
        }
    }

    if total_bases > 0.0 {
        let total_gc: f64 = stats
            .iter()
            .filter_map(|s| Some(s.get("gc_content")?.as_f64()? * s.get("total_bases")?.as_f64()?))
            .sum();
        merged.insert("gc_content".to_string(), Value::from(total_gc / total_bases));
    }
}

/// Merge read statistics objects (before/after filtering, or per-read sections).
///
/// `stage_name` controls the merge strategy:
/// - Names containing `"read"` (e.g. `"read1_before_filtering"`) are read-level sections:
///   `total_cycles` is copied from the first chunk; derived rates are not recalculated.
/// - All other names (e.g. `"before_filtering"`) are summary sections: rates and mean
///   lengths are recalculated from the merged totals.
#[must_use]
pub fn merge_read_stats(stats: &[Value], stage_name: &str) -> Value {
    let mut merged = Map::new();

    for field in &["total_reads", "total_bases", "q20_bases", "q30_bases"] {
        let sum: i64 = stats.iter().filter_map(|s| s.get(*field)?.as_i64()).sum();
        merged.insert((*field).to_string(), Value::from(sum));
    }

    if stage_name.contains("read") {
        if let Some(cycles) = stats.first().and_then(|s| s.get("total_cycles")) {
            merged.insert("total_cycles".to_string(), cycles.clone());
        }
    } else {
        apply_summary_rates(&mut merged, stats);
    }

    for (key, section) in &[
        ("quality_curves", "quality_curves"),
        ("content_curves", "content_curves"),
    ] {
        let maps: Vec<&Map<String, Value>> = stats
            .iter()
            .filter_map(|s| s.get(*section)?.as_object())
            .collect();
        if !maps.is_empty() {
            merged.insert((*key).to_string(), Value::Object(merge_curve_maps(&maps)));
        }
    }

    let kmer_counts: Vec<&Map<String, Value>> = stats
        .iter()
        .filter_map(|s| s.get("kmer_count")?.as_object())
        .collect();
    if !kmer_counts.is_empty() {
        merged.insert(
            "kmer_count".to_string(),
            Value::Object(merge_kmer_counts(&kmer_counts)),
        );
    }

    if stats.iter().any(|s| s.get("overrepresented_sequences").is_some()) {
        merged.insert(
            "overrepresented_sequences".to_string(),
            Value::Object(Map::new()),
        );
    }

    Value::Object(merged)
}

/// Merge adapter cutting objects by summing trimmed counts and merging per-adapter counts.
#[must_use]
pub fn merge_adapter_cutting(cuttings: &[Value]) -> Value {
    let mut merged = Map::new();

    let trimmed_reads: i64 = cuttings
        .iter()
        .filter_map(|ac| ac.get("adapter_trimmed_reads")?.as_i64())
        .sum();
    let trimmed_bases: i64 = cuttings
        .iter()
        .filter_map(|ac| ac.get("adapter_trimmed_bases")?.as_i64())
        .sum();
    merged.insert("adapter_trimmed_reads".to_string(), Value::from(trimmed_reads));
    merged.insert("adapter_trimmed_bases".to_string(), Value::from(trimmed_bases));

    // Adapter sequences are chunk-invariant; take from the first file. read2 fields
    // are absent for single-end data, so each is emitted only when actually present.
    for field in &["read1_adapter_sequence", "read2_adapter_sequence"] {
        if let Some(seq) = cuttings.first().and_then(|ac| ac.get(*field)) {
            merged.insert((*field).to_string(), seq.clone());
        }
    }

    for field in &["read1_adapter_counts", "read2_adapter_counts"] {
        if cuttings.iter().any(|ac| ac.get(*field).is_some()) {
            let counts: Vec<&Map<String, Value>> = cuttings
                .iter()
                .filter_map(|ac| ac.get(*field)?.as_object())
                .collect();
            merged.insert(
                (*field).to_string(),
                Value::Object(merge_adapter_counts(&counts)),
            );
        }
    }

    Value::Object(merged)
}

/// Merge a list of parsed fastp JSON objects into one, as if fastp had processed the whole dataset.
#[must_use]
pub fn merge_fastp_jsons(data_list: &[Value]) -> Value {
    if data_list.is_empty() {
        return Value::Object(Map::new());
    }

    let mut merged = Map::new();

    // Summary
    let mut summary = Map::new();
    if let Some(first_summary) = data_list[0].get("summary") {
        for field in &["fastp_version", "sequencing"] {
            if let Some(val) = first_summary.get(*field) {
                summary.insert((*field).to_string(), val.clone());
            }
        }
    }
    let before: Vec<Value> = data_list
        .iter()
        .filter_map(|d| d.get("summary")?.get("before_filtering").cloned())
        .collect();
    let after: Vec<Value> = data_list
        .iter()
        .filter_map(|d| d.get("summary")?.get("after_filtering").cloned())
        .collect();
    summary.insert(
        "before_filtering".to_string(),
        merge_read_stats(&before, "before_filtering"),
    );
    summary.insert(
        "after_filtering".to_string(),
        merge_read_stats(&after, "after_filtering"),
    );
    merged.insert("summary".to_string(), Value::Object(summary));

    // Top-level sections
    let filtering_results: Vec<Value> = data_list
        .iter()
        .filter_map(|d| d.get("filtering_result").cloned())
        .collect();
    merged.insert("filtering_result".to_string(), merge_filtering_result(&filtering_results));

    let dups: Vec<Value> = data_list
        .iter()
        .filter_map(|d| d.get("duplication").cloned())
        .collect();
    let summaries: Vec<Value> = data_list
        .iter()
        .filter_map(|d| d.get("summary").cloned())
        .collect();
    merged.insert("duplication".to_string(), merge_duplication(&dups, &summaries));

    // insert_size is paired-end only; absent for single-end data.
    if data_list[0].get("insert_size").is_some() {
        let insert_sizes: Vec<Value> = data_list
            .iter()
            .filter_map(|d| d.get("insert_size").cloned())
            .collect();
        merged.insert("insert_size".to_string(), merge_insert_size(&insert_sizes));
    }

    let adapter_cuttings: Vec<Value> = data_list
        .iter()
        .filter_map(|d| d.get("adapter_cutting").cloned())
        .collect();
    merged.insert(
        "adapter_cutting".to_string(),
        merge_adapter_cutting(&adapter_cuttings),
    );

    // Per-read sections (only present for paired-end data)
    for section in &[
        "read1_before_filtering",
        "read1_after_filtering",
        "read2_before_filtering",
        "read2_after_filtering",
    ] {
        if data_list[0].get(*section).is_some() {
            let values: Vec<Value> = data_list
                .iter()
                .filter_map(|d| d.get(*section).cloned())
                .collect();
            merged.insert((*section).to_string(), merge_read_stats(&values, section));
        }
    }

    Value::Object(merged)
}
