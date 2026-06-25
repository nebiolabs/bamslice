//! Merge multiple fastp JSON output files into a single aggregated result,
//! as if fastp had been run on the entire dataset at once.
//!
//! The merge fails loudly: a required field that is missing or has the wrong
//! type produces an error rather than being silently coerced to zero, so a
//! malformed input never yields a plausible-but-wrong merged result. Genuinely
//! optional content (single-end data lacks `read2_*` sections and `insert_size`,
//! for example) is skipped when absent but still validated when present.

use anyhow::{Context, Result};
use serde_json::{Map, Value};

/// Extract a required field as `i64`, erroring if it is absent or not an integer.
fn required_i64(obj: &Value, field: &str) -> Result<i64> {
    let value = obj
        .get(field)
        .with_context(|| format!("missing required field '{field}'"))?;
    value
        .as_i64()
        .with_context(|| format!("field '{field}' is not an integer (found {value})"))
}

/// Extract a required field as `f64`, erroring if it is absent or not a number.
fn required_f64(obj: &Value, field: &str) -> Result<f64> {
    let value = obj
        .get(field)
        .with_context(|| format!("missing required field '{field}'"))?;
    value
        .as_f64()
        .with_context(|| format!("field '{field}' is not a number (found {value})"))
}

/// Coerce an array element to `f64`, erroring if it is not a number.
fn element_as_f64(value: &Value) -> Result<f64> {
    value
        .as_f64()
        .with_context(|| format!("array contains a non-numeric value ({value})"))
}

/// Sum arrays element-wise, padding shorter arrays with zero.
///
/// Preserves integer values when all inputs are integers; uses `f64` otherwise.
///
/// # Errors
///
/// Returns an error if any present element is not a number.
pub fn sum_arrays(arrays: &[&[Value]]) -> Result<Vec<Value>> {
    if arrays.is_empty() {
        return Ok(Vec::new());
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
        Ok(result.into_iter().map(Value::from).collect())
    } else {
        let mut result = vec![0_f64; max_len];
        for arr in arrays {
            for (i, val) in arr.iter().enumerate() {
                result[i] += element_as_f64(val)?;
            }
        }
        Ok(result.into_iter().map(Value::from).collect())
    }
}

/// Average arrays element-wise, weighting each array `i` by `weights[i]`.
///
/// fastp reports each per-cycle curve value as a mean over a chunk's reads, so the
/// whole-dataset mean is the read-count-weighted average of the per-chunk means.
/// Each position accumulates only the weights of arrays that actually reach it, so
/// variable-length inputs are handled without padding bias. A position with zero
/// total weight yields `0.0`.
///
/// # Errors
///
/// Returns an error if any present element is not a number.
pub fn weighted_average_arrays(arrays: &[&[Value]], weights: &[f64]) -> Result<Vec<Value>> {
    if arrays.is_empty() {
        return Ok(Vec::new());
    }
    let max_len = arrays.iter().map(|a| a.len()).max().unwrap_or(0);
    let mut sums = vec![0_f64; max_len];
    let mut weight_totals = vec![0_f64; max_len];
    for (arr, &weight) in arrays.iter().zip(weights.iter()) {
        for (i, val) in arr.iter().enumerate() {
            sums[i] = element_as_f64(val)?.mul_add(weight, sums[i]);
            weight_totals[i] += weight;
        }
    }
    Ok(sums
        .iter()
        .zip(weight_totals.iter())
        .map(|(&sum, &total)| {
            if total > 0.0 {
                Value::from(sum / total)
            } else {
                Value::from(0.0)
            }
        })
        .collect())
}

/// Merge string-keyed integer count maps by summing per-key counts.
///
/// `item_label` is used in error messages to identify the key type (e.g. "kmer", "adapter").
/// Errors if any count is not an integer.
fn merge_int_counts(dicts: &[&Map<String, Value>], item_label: &str) -> Result<Map<String, Value>> {
    let mut merged: Map<String, Value> = Map::new();
    for dict in dicts {
        for (key, count) in *dict {
            let count = count.as_i64().with_context(|| {
                format!("{item_label} '{key}' count is not an integer (found {count})")
            })?;
            let existing = merged.get(key).and_then(Value::as_i64).unwrap_or(0);
            merged.insert(key.clone(), Value::from(existing.saturating_add(count)));
        }
    }
    Ok(merged)
}

fn merge_kmer_counts(dicts: &[&Map<String, Value>]) -> Result<Map<String, Value>> {
    merge_int_counts(dicts, "kmer")
}

/// Merge quality or content curve objects by averaging each key's array across chunks,
/// weighting every chunk by its read count (`weights[i]`, index-aligned with `curves`).
///
/// A chunk's curve value at a cycle is a mean over that chunk's reads, so the
/// dataset-wide mean is the read-count-weighted average of the per-chunk means — not a
/// flat average, which would over-weight small chunks. `weights` therefore carries each
/// chunk's section `total_reads`. See [`weighted_average_arrays`] for the per-cycle detail.
///
/// Errors if a curve value present in a chunk is not an array of numbers.
fn merge_curve_maps(curves: &[&Map<String, Value>], weights: &[f64]) -> Result<Map<String, Value>> {
    if curves.is_empty() {
        return Ok(Map::new());
    }
    let mut seen = std::collections::HashSet::new();
    let all_keys: Vec<&String> = curves
        .iter()
        .flat_map(|c| c.keys())
        .filter(|k| seen.insert(k.as_str()))
        .collect();

    let mut merged = Map::new();
    for key in all_keys {
        let mut arrays: Vec<&[Value]> = Vec::new();
        let mut arr_weights: Vec<f64> = Vec::new();
        for (curve, &weight) in curves.iter().zip(weights.iter()) {
            if let Some(value) = curve.get(key) {
                let array = value
                    .as_array()
                    .with_context(|| format!("curve '{key}' is not an array"))?;
                arrays.push(array.as_slice());
                arr_weights.push(weight);
            }
        }
        let averaged = weighted_average_arrays(&arrays, &arr_weights)
            .with_context(|| format!("averaging curve '{key}'"))?;
        merged.insert(key.clone(), Value::Array(averaged));
    }
    Ok(merged)
}

fn merge_adapter_counts(dicts: &[&Map<String, Value>]) -> Result<Map<String, Value>> {
    merge_int_counts(dicts, "adapter")
}

/// Merge filtering result objects by summing all count fields.
///
/// Errors if any required count field is missing or not an integer.
fn merge_filtering_result(results: &[Value]) -> Result<Value> {
    const FIELDS: &[&str] = &[
        "passed_filter_reads",
        "low_quality_reads",
        "too_many_N_reads",
        "too_short_reads",
        "too_long_reads",
    ];
    let mut merged = Map::new();
    for field in FIELDS {
        let mut sum = 0_i64;
        for (i, result) in results.iter().enumerate() {
            sum += required_i64(result, field).with_context(|| format!("filtering_result[{i}]"))?;
        }
        merged.insert((*field).to_string(), Value::from(sum));
    }
    Ok(Value::Object(merged))
}

/// Merge insert size objects by summing histograms and recomputing the peak.
///
/// Note: the summing here is exact, but read pairs split across a chunk boundary
/// lose their mate and get counted as `unknown` by fastp, so the merged result drifts
/// slightly from a whole-dataset run. That is a slicing artifact, not a merge bug — see
/// `docs/merge-fidelity.md`.
///
/// Errors if a histogram or the `unknown` count is missing or malformed.
fn merge_insert_size(sizes: &[Value]) -> Result<Value> {
    let mut histograms: Vec<&[Value]> = Vec::new();
    for (i, size) in sizes.iter().enumerate() {
        let histogram = size
            .get("histogram")
            .with_context(|| format!("insert_size[{i}]: missing 'histogram'"))?
            .as_array()
            .with_context(|| format!("insert_size[{i}]: 'histogram' is not an array"))?;
        histograms.push(histogram.as_slice());
    }
    let merged_histogram = sum_arrays(&histograms).context("insert_size histogram")?;

    let mut unknown = 0_i64;
    for (i, size) in sizes.iter().enumerate() {
        unknown += required_i64(size, "unknown").with_context(|| format!("insert_size[{i}]"))?;
    }

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
    Ok(Value::Object(merged))
}

/// Compute a weighted-average duplication rate, weighting each chunk by its total read count.
///
/// This is an approximation: duplicates that fall in different chunks are invisible to
/// any per-chunk statistic, so the merged rate undercounts cross-chunk duplicates and
/// cannot match a whole-dataset run. See `docs/merge-fidelity.md`.
///
/// Errors if a chunk's duplication rate or its `summary.before_filtering.total_reads`
/// (the weight) is missing or malformed.
fn merge_duplication(dups: &[Value], summaries: &[Value]) -> Result<Value> {
    let mut weighted_sum = 0.0_f64;
    let mut total_reads = 0.0_f64;
    for (i, (dup, summary)) in dups.iter().zip(summaries.iter()).enumerate() {
        let reads = summary
            .get("before_filtering")
            .and_then(|b| b.get("total_reads"))
            .with_context(|| {
                format!("duplication[{i}]: missing summary.before_filtering.total_reads")
            })?
            .as_f64()
            .with_context(|| format!("duplication[{i}]: total_reads is not a number"))?;
        let rate = required_f64(dup, "rate").with_context(|| format!("duplication[{i}]"))?;
        weighted_sum = rate.mul_add(reads, weighted_sum);
        total_reads += reads;
    }
    let rate = if total_reads > 0.0 {
        weighted_sum / total_reads
    } else {
        0.0
    };
    let mut merged = Map::new();
    merged.insert("rate".to_string(), Value::from(rate));
    Ok(Value::Object(merged))
}

/// Recalculate derived rate and mean-length fields for a summary section after totals are merged.
/// Only called for `before_filtering` / `after_filtering` sections, not per-read sections.
#[allow(clippy::cast_possible_truncation)]
fn apply_summary_rates(merged: &mut Map<String, Value>, stats: &[Value]) -> Result<()> {
    // total_bases / total_reads were just inserted by the caller as integers.
    let total_bases = merged
        .get("total_bases")
        .and_then(Value::as_f64)
        .unwrap_or(0.0);
    let total_reads = merged
        .get("total_reads")
        .and_then(Value::as_f64)
        .unwrap_or(0.0);

    if total_bases > 0.0 {
        let q20 = merged
            .get("q20_bases")
            .and_then(Value::as_f64)
            .unwrap_or(0.0);
        let q30 = merged
            .get("q30_bases")
            .and_then(Value::as_f64)
            .unwrap_or(0.0);
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
            if stats.iter().all(|s| s.get(*field).is_none()) {
                continue;
            }
            let mut weighted_sum = 0.0_f64;
            for (i, stat) in stats.iter().enumerate() {
                let mean = required_f64(stat, field).with_context(|| format!("summary[{i}]"))?;
                let reads =
                    required_f64(stat, "total_reads").with_context(|| format!("summary[{i}]"))?;
                weighted_sum = mean.mul_add(reads, weighted_sum);
            }
            // Truncation is intentional: fastp stores mean lengths as integers
            merged.insert(
                (*field).to_string(),
                Value::from((weighted_sum / total_reads) as i64),
            );
        }
    }

    if total_bases > 0.0 && stats.iter().any(|s| s.get("gc_content").is_some()) {
        let mut total_gc = 0.0_f64;
        for (i, stat) in stats.iter().enumerate() {
            let gc = required_f64(stat, "gc_content").with_context(|| format!("summary[{i}]"))?;
            let bases =
                required_f64(stat, "total_bases").with_context(|| format!("summary[{i}]"))?;
            total_gc = gc.mul_add(bases, total_gc);
        }
        merged.insert(
            "gc_content".to_string(),
            Value::from(total_gc / total_bases),
        );
    }

    Ok(())
}

/// Distinguishes summary sections (`before_filtering`, `after_filtering`) from per-read
/// sections (`read1_before_filtering`, etc.) when calling [`merge_read_stats`].
#[derive(Clone, Copy)]
pub enum ReadStatsKind {
    /// Rates and mean lengths are recomputed from merged totals.
    Summary,
    /// `total_cycles` is copied from the first chunk; rates are not recalculated.
    PerRead,
}

/// Merge read statistics objects (before/after filtering, or per-read sections).
///
/// # Errors
///
/// Returns an error if a required count field is missing or not an integer, or if an
/// optional sub-section (curves, kmer counts) is present but malformed.
pub fn merge_read_stats(stats: &[Value], stage_name: &str, kind: ReadStatsKind) -> Result<Value> {
    let mut merged = Map::new();

    for field in &["total_reads", "total_bases", "q20_bases", "q30_bases"] {
        let mut sum = 0_i64;
        for (i, stat) in stats.iter().enumerate() {
            sum += required_i64(stat, field).with_context(|| format!("{stage_name}[{i}]"))?;
        }
        merged.insert((*field).to_string(), Value::from(sum));
    }

    match kind {
        ReadStatsKind::PerRead => {
            if let Some(cycles) = stats.first().and_then(|s| s.get("total_cycles")) {
                merged.insert("total_cycles".to_string(), cycles.clone());
            }
        }
        ReadStatsKind::Summary => {
            apply_summary_rates(&mut merged, stats)?;
        }
    }

    for section in &["quality_curves", "content_curves"] {
        // Curves are read-count-weighted means, so each chunk's curve is paired with
        // that chunk's `total_reads` (the per-cycle weight; exact for raw, full-length
        // reads, a close proxy once trimming makes lengths vary).
        let mut maps: Vec<&Map<String, Value>> = Vec::new();
        let mut weights: Vec<f64> = Vec::new();
        for stat in stats {
            if let Some(obj) = stat.get(*section).and_then(Value::as_object) {
                maps.push(obj);
                weights.push(
                    stat.get("total_reads")
                        .and_then(Value::as_f64)
                        .unwrap_or(0.0),
                );
            }
        }
        if !maps.is_empty() {
            let curves = merge_curve_maps(&maps, &weights)
                .with_context(|| format!("{stage_name} {section}"))?;
            merged.insert((*section).to_string(), Value::Object(curves));
        }
    }

    let kmer_counts: Vec<&Map<String, Value>> = stats
        .iter()
        .filter_map(|s| s.get("kmer_count")?.as_object())
        .collect();
    if !kmer_counts.is_empty() {
        let merged_kmers =
            merge_kmer_counts(&kmer_counts).with_context(|| format!("{stage_name} kmer_count"))?;
        merged.insert("kmer_count".to_string(), Value::Object(merged_kmers));
    }

    // overrepresented_sequences cannot be meaningfully merged: each chunk prunes its
    // own top-N sequences independently, so the per-chunk sets are not additive and
    // the discarded tail can't be recovered. Emit an empty object (fastp's own output
    // when nothing is found) to preserve the key's presence in the merged JSON.
    if stats
        .iter()
        .any(|s| s.get("overrepresented_sequences").is_some())
    {
        merged.insert(
            "overrepresented_sequences".to_string(),
            Value::Object(Map::new()),
        );
    }

    Ok(Value::Object(merged))
}

/// Merge adapter cutting objects by summing trimmed counts and merging per-adapter counts.
///
/// Errors if the trimmed-read/base counts are missing or not integers. read2 fields are
/// absent for single-end data, so each is emitted only when present (and validated if so).
fn merge_adapter_cutting(cuttings: &[Value]) -> Result<Value> {
    let mut merged = Map::new();

    let mut trimmed_reads = 0_i64;
    let mut trimmed_bases = 0_i64;
    for (i, cutting) in cuttings.iter().enumerate() {
        trimmed_reads += required_i64(cutting, "adapter_trimmed_reads")
            .with_context(|| format!("adapter_cutting[{i}]"))?;
        trimmed_bases += required_i64(cutting, "adapter_trimmed_bases")
            .with_context(|| format!("adapter_cutting[{i}]"))?;
    }
    merged.insert(
        "adapter_trimmed_reads".to_string(),
        Value::from(trimmed_reads),
    );
    merged.insert(
        "adapter_trimmed_bases".to_string(),
        Value::from(trimmed_bases),
    );

    // Adapter sequences are chunk-invariant; take from the first file. read2 fields
    // are absent for single-end data, so each is emitted only when actually present.
    for field in &["read1_adapter_sequence", "read2_adapter_sequence"] {
        if let Some(seq) = cuttings.first().and_then(|ac| ac.get(*field)) {
            merged.insert((*field).to_string(), seq.clone());
        }
    }

    for field in &["read1_adapter_counts", "read2_adapter_counts"] {
        if cuttings.iter().any(|ac| ac.get(*field).is_some()) {
            let mut counts: Vec<&Map<String, Value>> = Vec::new();
            for (i, cutting) in cuttings.iter().enumerate() {
                if let Some(value) = cutting.get(*field) {
                    let object = value.as_object().with_context(|| {
                        format!("adapter_cutting[{i}]: '{field}' is not an object")
                    })?;
                    counts.push(object);
                }
            }
            let merged_counts = merge_adapter_counts(&counts)
                .with_context(|| format!("adapter_cutting {field}"))?;
            merged.insert((*field).to_string(), Value::Object(merged_counts));
        }
    }

    Ok(Value::Object(merged))
}

/// Extract a required nested section from every input file, erroring (with the file
/// index and the dotted path) if any file is missing it.
fn require_section(data_list: &[Value], keys: &[&str]) -> Result<Vec<Value>> {
    let mut out = Vec::with_capacity(data_list.len());
    for (i, data) in data_list.iter().enumerate() {
        let mut current = data;
        for key in keys {
            current = current.get(*key).with_context(|| {
                format!("file[{i}]: missing required section '{}'", keys.join("."))
            })?;
        }
        out.push(current.clone());
    }
    Ok(out)
}

/// Merge a list of parsed fastp JSON objects into one, as if fastp had processed the whole dataset.
///
/// # Errors
///
/// Returns an error if any input file is missing a required section
/// (`summary.before_filtering`, `summary.after_filtering`, `filtering_result`,
/// `duplication`, `adapter_cutting`), or if any required numeric field within a
/// merged section is missing or has a non-numeric value.
pub fn merge_jsons(data_list: &[Value]) -> Result<Value> {
    if data_list.is_empty() {
        return Ok(Value::Object(Map::new()));
    }

    let mut merged = Map::new();

    // Summary: copy chunk-invariant metadata from the first file, merge the rest.
    let mut summary = Map::new();
    if let Some(first_summary) = data_list[0].get("summary") {
        for field in &["fastp_version", "sequencing"] {
            if let Some(val) = first_summary.get(*field) {
                summary.insert((*field).to_string(), val.clone());
            }
        }
    }
    let before = require_section(data_list, &["summary", "before_filtering"])?;
    let after = require_section(data_list, &["summary", "after_filtering"])?;
    summary.insert(
        "before_filtering".to_string(),
        merge_read_stats(&before, "before_filtering", ReadStatsKind::Summary)?,
    );
    summary.insert(
        "after_filtering".to_string(),
        merge_read_stats(&after, "after_filtering", ReadStatsKind::Summary)?,
    );
    merged.insert("summary".to_string(), Value::Object(summary));

    // Top-level sections present for both single- and paired-end data.
    let filtering_results = require_section(data_list, &["filtering_result"])?;
    merged.insert(
        "filtering_result".to_string(),
        merge_filtering_result(&filtering_results)?,
    );

    let dups = require_section(data_list, &["duplication"])?;
    let summaries = require_section(data_list, &["summary"])?;
    merged.insert(
        "duplication".to_string(),
        merge_duplication(&dups, &summaries)?,
    );

    // insert_size is paired-end only; absent for single-end data. Emitted before
    // adapter_cutting to match fastp's own top-level key ordering.
    if data_list.iter().any(|d| d.get("insert_size").is_some()) {
        let insert_sizes = require_section(data_list, &["insert_size"])?;
        merged.insert("insert_size".to_string(), merge_insert_size(&insert_sizes)?);
    }

    let adapter_cuttings = require_section(data_list, &["adapter_cutting"])?;
    merged.insert(
        "adapter_cutting".to_string(),
        merge_adapter_cutting(&adapter_cuttings)?,
    );

    // Per-read sections (only present for paired-end data, except read1 for single-end).
    // Order matches fastp's own output: both `before` sections, then both `after`.
    for section in &[
        "read1_before_filtering",
        "read2_before_filtering",
        "read1_after_filtering",
        "read2_after_filtering",
    ] {
        if data_list.iter().any(|d| d.get(*section).is_some()) {
            let values = require_section(data_list, &[section])?;
            merged.insert(
                (*section).to_string(),
                merge_read_stats(&values, section, ReadStatsKind::PerRead)?,
            );
        }
    }

    Ok(Value::Object(merged))
}

/// Strip trailing zeros (and a now-orphaned decimal point) from a fixed-notation
/// number, matching `%g`'s suppression of insignificant trailing digits.
fn strip_trailing_zeros(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

/// Format an `f64` the way fastp's C++ output does.
///
/// fastp writes its JSON with a `std::ostream` left at the default precision, which
/// is `printf`'s `%g` with 6 significant figures. Reproducing its bytes therefore
/// means: round to 6 significant digits; use fixed notation when the decimal
/// exponent `X` satisfies `-4 <= X < 6` and scientific (`d.ddddde±XX`, exponent sign
/// always present and zero-padded to at least two digits) otherwise; strip trailing
/// zeros and any orphaned decimal point. Zero is written as a bare `0`.
fn fastp_float_repr(f: f64) -> String {
    if f == 0.0 {
        return "0".to_string();
    }
    if f.is_nan() {
        return "nan".to_string();
    }
    if f.is_infinite() {
        return if f < 0.0 { "-inf" } else { "inf" }.to_string();
    }

    // `{:.5e}` yields exactly 6 significant digits in scientific form
    // (e.g. "4.12386e-1"), performing the round-to-6-sig-figs for us.
    let sci = format!("{f:.5e}");
    let (mantissa, exp_str) = sci
        .split_once('e')
        .expect("scientific format always contains 'e'");
    let exp: i32 = exp_str.parse().expect("scientific exponent is an integer");

    let negative = mantissa.starts_with('-');
    let mant = mantissa.trim_start_matches('-');
    // The six significant digits, decimal point removed.
    let digits: String = mant.chars().filter(|&c| c != '.').collect();

    let mut out = String::new();
    if negative {
        out.push('-');
    }

    // %g uses fixed notation when -4 <= exp < precision (6), scientific otherwise.
    if (-4..6).contains(&exp) {
        // Position of the decimal point relative to the start of `digits`.
        let decpt = exp + 1;
        let body = if decpt <= 0 {
            let mut s = String::from("0.");
            for _ in 0..-decpt {
                s.push('0');
            }
            s.push_str(&digits);
            s
        } else {
            let split = usize::try_from(decpt).unwrap_or(0);
            if split >= digits.len() {
                let mut s = digits.clone();
                for _ in 0..(split - digits.len()) {
                    s.push('0');
                }
                s
            } else {
                format!("{}.{}", &digits[..split], &digits[split..])
            }
        };
        out.push_str(strip_trailing_zeros(&body));
    } else {
        // Scientific: d.ddddd, with trailing zeros stripped from the fraction.
        let rest = digits[1..].trim_end_matches('0');
        out.push_str(&digits[..1]);
        if !rest.is_empty() {
            out.push('.');
            out.push_str(rest);
        }
        out.push('e');
        out.push(if exp < 0 { '-' } else { '+' });
        let magnitude = exp.abs();
        if magnitude < 10 {
            out.push('0');
        }
        out.push_str(&magnitude.to_string());
    }

    out
}

/// Format a JSON number, preserving integers verbatim and rendering floats as fastp does.
fn format_number(n: &serde_json::Number) -> String {
    if n.is_i64() || n.is_u64() {
        n.to_string()
    } else if let Some(f) = n.as_f64() {
        fastp_float_repr(f)
    } else {
        n.to_string()
    }
}

/// Append one indentation level (`depth` tabs) to `out`.
fn push_indent(out: &mut String, depth: usize) {
    for _ in 0..depth {
        out.push('\t');
    }
}

/// Write an array compacted onto a single line: `[v,v,v]`, no interior whitespace.
/// fastp arrays are numeric; any non-number element is written compactly as well.
fn write_array(arr: &[Value], out: &mut String) -> Result<()> {
    out.push('[');
    for (i, v) in arr.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        match v {
            Value::Number(n) => out.push_str(&format_number(n)),
            other => write_value(other, 0, out)?,
        }
    }
    out.push(']');
    Ok(())
}

/// Write a `kmer_count` object as a matrix of 16 `"KMER": count` entries per line,
/// matching fastp's layout. `depth` is the indentation of the enclosing braces.
fn write_kmer_count(map: &Map<String, Value>, depth: usize, out: &mut String) -> Result<()> {
    if map.is_empty() {
        out.push_str("{}");
        return Ok(());
    }
    let entries: Vec<String> = map
        .iter()
        .map(|(k, v)| {
            Ok(format!(
                "{}: {}",
                serde_json::to_string(k)?,
                format_number_value(v)
            ))
        })
        .collect::<Result<_>>()?;

    out.push_str("{\n");
    let mut i = 0;
    while i < entries.len() {
        let end = (i + 16).min(entries.len());
        push_indent(out, depth + 1);
        out.push_str(&entries[i..end].join(",\t\t\t"));
        if end < entries.len() {
            out.push(',');
        }
        out.push('\n');
        i += 16;
    }
    push_indent(out, depth);
    out.push('}');
    Ok(())
}

/// Format a `Value` known to be a number (kmer counts are always integers).
fn format_number_value(v: &Value) -> String {
    match v {
        Value::Number(n) => format_number(n),
        other => other.to_string(),
    }
}

/// Write an object across multiple lines with tab indentation. Array-valued members
/// get no space after the colon (`"key":[...]`) to match fastp; everything else does.
fn write_object(map: &Map<String, Value>, depth: usize, out: &mut String) -> Result<()> {
    if map.is_empty() {
        out.push_str("{}");
        return Ok(());
    }
    out.push_str("{\n");
    let n = map.len();
    for (i, (k, v)) in map.iter().enumerate() {
        push_indent(out, depth + 1);
        out.push_str(&serde_json::to_string(k)?);
        match v {
            Value::Array(arr) => {
                out.push(':');
                write_array(arr, out)?;
            }
            Value::Object(obj) if k == "kmer_count" => {
                out.push_str(": ");
                write_kmer_count(obj, depth + 1, out)?;
            }
            _ => {
                out.push_str(": ");
                write_value(v, depth + 1, out)?;
            }
        }
        if i + 1 < n {
            out.push(',');
        }
        out.push('\n');
    }
    push_indent(out, depth);
    out.push('}');
    Ok(())
}

/// Recursively write a `Value` in fastp's JSON layout.
fn write_value(value: &Value, depth: usize, out: &mut String) -> Result<()> {
    match value {
        Value::Object(map) => write_object(map, depth, out),
        Value::Array(arr) => write_array(arr, out),
        Value::String(_) => {
            out.push_str(&serde_json::to_string(value)?);
            Ok(())
        }
        Value::Number(n) => {
            out.push_str(&format_number(n));
            Ok(())
        }
        Value::Bool(b) => {
            out.push_str(if *b { "true" } else { "false" });
            Ok(())
        }
        Value::Null => {
            out.push_str("null");
            Ok(())
        }
    }
}

/// Serialize a merged fastp value to match fastp's own JSON layout byte-for-byte.
///
/// Objects are pretty-printed with tab indentation, arrays are compacted onto a
/// single line with no interior whitespace and no space after the colon
/// (`"histogram":[1,2,3]`), `kmer_count` maps are laid out as a matrix of 16 entries
/// per line, and floats are formatted as fastp's C++ output does (6 significant
/// figures; see [`fastp_float_repr`]). This reproduces fastp's own layout so the
/// result diffs cleanly against an unsplit run.
///
/// The returned string has no trailing newline; callers append one to match fastp.
///
/// # Errors
///
/// Returns an error if a key or string value cannot be serialized.
pub fn format_json(value: &Value) -> Result<String> {
    let mut out = String::new();
    write_value(value, 0, &mut out)?;
    Ok(out)
}
