use bamslice::fastp_merge::{
    average_arrays, merge_adapter_cutting, merge_adapter_counts, merge_curve_maps,
    merge_duplication, merge_fastp_jsons, merge_filtering_result, merge_insert_size,
    merge_kmer_counts, merge_read_stats, sum_arrays,
};
use serde_json::{json, Value};
use std::fs;

fn load_fixture(name: &str) -> Value {
    let path = format!("tests/fixtures/{name}");
    let content = fs::read_to_string(&path).unwrap_or_else(|e| panic!("Cannot read {path}: {e}"));
    serde_json::from_str(&content)
        .unwrap_or_else(|e| panic!("Cannot parse JSON from {path}: {e}"))
}

fn assert_approx_eq(actual: f64, expected: f64, label: &str) {
    let diff = (actual - expected).abs();
    assert!(
        diff < 1e-9,
        "{label}: expected {expected}, got {actual} (diff={diff})"
    );
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[test]
fn test_sum_arrays_integers() {
    let a = vec![json!(1), json!(2), json!(3)];
    let b = vec![json!(4), json!(5), json!(6)];
    let result = sum_arrays(&[a.as_slice(), b.as_slice()]);
    assert_eq!(result, vec![json!(5), json!(7), json!(9)]);
}

#[test]
fn test_sum_arrays_variable_length() {
    let a = vec![json!(1), json!(2)];
    let b = vec![json!(3), json!(4), json!(5)];
    let result = sum_arrays(&[a.as_slice(), b.as_slice()]);
    assert_eq!(result, vec![json!(4), json!(6), json!(5)]);
}

#[test]
fn test_sum_arrays_empty() {
    let result = sum_arrays(&[]);
    assert!(result.is_empty());
}

#[test]
fn test_average_arrays() {
    let a = vec![json!(10.0), json!(20.0)];
    let b = vec![json!(20.0), json!(40.0)];
    let result = average_arrays(&[a.as_slice(), b.as_slice()]);
    assert_approx_eq(result[0].as_f64().unwrap(), 15.0, "index 0");
    assert_approx_eq(result[1].as_f64().unwrap(), 30.0, "index 1");
}

#[test]
fn test_average_arrays_variable_length() {
    let a = vec![json!(10.0)];
    let b = vec![json!(20.0), json!(30.0)];
    let result = average_arrays(&[a.as_slice(), b.as_slice()]);
    assert_approx_eq(result[0].as_f64().unwrap(), 15.0, "index 0");
    assert_approx_eq(result[1].as_f64().unwrap(), 30.0, "index 1 (only one source)");
}

#[test]
fn test_merge_kmer_counts() {
    let mut d1 = serde_json::Map::new();
    d1.insert("AAAAA".to_string(), json!(100));
    d1.insert("AAAAC".to_string(), json!(50));
    let mut d2 = serde_json::Map::new();
    d2.insert("AAAAA".to_string(), json!(200));
    d2.insert("AAAAT".to_string(), json!(75));
    let result = merge_kmer_counts(&[&d1, &d2]);
    assert_eq!(result["AAAAA"].as_i64(), Some(300));
    assert_eq!(result["AAAAC"].as_i64(), Some(50));
    assert_eq!(result["AAAAT"].as_i64(), Some(75));
}

#[test]
fn test_merge_curve_maps() {
    let mut c1 = serde_json::Map::new();
    c1.insert("mean".to_string(), json!([30.0, 35.0, 40.0]));
    let mut c2 = serde_json::Map::new();
    c2.insert("mean".to_string(), json!([40.0, 45.0, 50.0]));
    let result = merge_curve_maps(&[&c1, &c2]);
    let mean = result["mean"].as_array().unwrap();
    assert_approx_eq(mean[0].as_f64().unwrap(), 35.0, "mean[0]");
    assert_approx_eq(mean[1].as_f64().unwrap(), 40.0, "mean[1]");
    assert_approx_eq(mean[2].as_f64().unwrap(), 45.0, "mean[2]");
}

#[test]
fn test_merge_adapter_counts() {
    let mut d1 = serde_json::Map::new();
    d1.insert("ADAPTERX".to_string(), json!(100));
    let mut d2 = serde_json::Map::new();
    d2.insert("ADAPTERX".to_string(), json!(200));
    d2.insert("ADAPTERY".to_string(), json!(50));
    let result = merge_adapter_counts(&[&d1, &d2]);
    assert_eq!(result["ADAPTERX"].as_i64(), Some(300));
    assert_eq!(result["ADAPTERY"].as_i64(), Some(50));
}

#[test]
fn test_merge_filtering_result() {
    let r1 = json!({
        "passed_filter_reads": 950, "low_quality_reads": 30,
        "too_many_N_reads": 10, "too_short_reads": 10, "too_long_reads": 0
    });
    let r2 = json!({
        "passed_filter_reads": 1900, "low_quality_reads": 60,
        "too_many_N_reads": 20, "too_short_reads": 20, "too_long_reads": 0
    });
    let result = merge_filtering_result(&[r1, r2]);
    assert_eq!(result["passed_filter_reads"].as_i64(), Some(2850));
    assert_eq!(result["low_quality_reads"].as_i64(), Some(90));
    assert_eq!(result["too_many_N_reads"].as_i64(), Some(30));
    assert_eq!(result["too_short_reads"].as_i64(), Some(30));
    assert_eq!(result["too_long_reads"].as_i64(), Some(0));
}

#[test]
fn test_merge_insert_size_sums_histogram_and_recomputes_peak() {
    let s1 = json!({"peak": 3, "unknown": 50, "histogram": [0, 5, 10, 20, 50, 100, 150, 200, 150, 100, 50]});
    let s2 = json!({"peak": 3, "unknown": 100, "histogram": [0, 10, 20, 40, 100, 200, 300, 400, 300, 200, 100]});
    let result = merge_insert_size(&[s1, s2]);
    // Merged histogram: [0,15,30,60,150,300,450,600,450,300,150]
    let hist = result["histogram"].as_array().unwrap();
    assert_eq!(hist[0].as_i64(), Some(0));
    assert_eq!(hist[7].as_i64(), Some(600)); // max → peak
    assert_eq!(hist[10].as_i64(), Some(150));
    assert_eq!(result["peak"].as_u64(), Some(7));
    assert_eq!(result["unknown"].as_i64(), Some(150));
}

#[test]
fn test_merge_duplication_weighted_average() {
    // chunk1: rate=0.1, reads=1000; chunk2: rate=0.2, reads=2000
    // expected: (0.1*1000 + 0.2*2000) / 3000 = 500/3000 ≈ 0.16667
    let d1 = json!({"rate": 0.1});
    let d2 = json!({"rate": 0.2});
    let s1 = json!({"before_filtering": {"total_reads": 1000}});
    let s2 = json!({"before_filtering": {"total_reads": 2000}});
    let result = merge_duplication(&[d1, d2], &[s1, s2]);
    assert_approx_eq(result["rate"].as_f64().unwrap(), 500.0 / 3000.0, "rate");
}

#[test]
fn test_merge_read_stats_summary_recalculates_rates() {
    let s1 = json!({
        "total_reads": 1000, "total_bases": 150_000,
        "q20_bases": 140_000, "q30_bases": 120_000,
        "q20_rate": 0.933, "q30_rate": 0.8,
        "read1_mean_length": 150, "read2_mean_length": 150,
        "gc_content": 0.42
    });
    let s2 = json!({
        "total_reads": 2000, "total_bases": 300_000,
        "q20_bases": 280_000, "q30_bases": 240_000,
        "q20_rate": 0.933, "q30_rate": 0.8,
        "read1_mean_length": 150, "read2_mean_length": 150,
        "gc_content": 0.40
    });
    let result = merge_read_stats(&[s1, s2], "before_filtering");

    assert_eq!(result["total_reads"].as_i64(), Some(3000));
    assert_eq!(result["total_bases"].as_i64(), Some(450_000));
    assert_eq!(result["q20_bases"].as_i64(), Some(420_000));
    // q20_rate = 420000/450000
    assert_approx_eq(result["q20_rate"].as_f64().unwrap(), 420_000.0 / 450_000.0, "q20_rate");
    assert_approx_eq(result["q30_rate"].as_f64().unwrap(), 360_000.0 / 450_000.0, "q30_rate");
    // gc_content = (0.42*150000 + 0.40*300000) / 450000
    assert_approx_eq(result["gc_content"].as_f64().unwrap(), 183_000.0 / 450_000.0, "gc_content");
    assert_eq!(result["read1_mean_length"].as_i64(), Some(150));
    assert_eq!(result["read2_mean_length"].as_i64(), Some(150));
    // total_cycles must NOT be present in summary sections
    assert!(result.get("total_cycles").is_none());
}

#[test]
fn test_merge_read_stats_read_section_preserves_cycles() {
    let s1 = json!({"total_reads": 500, "total_bases": 75000, "q20_bases": 70000, "q30_bases": 60000, "total_cycles": 150});
    let s2 = json!({"total_reads": 500, "total_bases": 75000, "q20_bases": 70000, "q30_bases": 60000, "total_cycles": 150});
    let result = merge_read_stats(&[s1, s2], "read1_before_filtering");

    assert_eq!(result["total_reads"].as_i64(), Some(1000));
    assert_eq!(result["total_cycles"].as_i64(), Some(150));
    // Rates must NOT be recalculated for read-level sections
    assert!(result.get("q20_rate").is_none());
    assert!(result.get("gc_content").is_none());
}

#[test]
fn test_merge_adapter_cutting() {
    let a1 = json!({
        "adapter_trimmed_reads": 200, "adapter_trimmed_bases": 5000,
        "read1_adapter_sequence": "AGATCGGAAGAGC",
        "read2_adapter_sequence": "AGATCGGAAGAGC",
        "read1_adapter_counts": {"AGATCGGAAGAGC": 150},
        "read2_adapter_counts": {"AGATCGGAAGAGC": 140}
    });
    let a2 = json!({
        "adapter_trimmed_reads": 400, "adapter_trimmed_bases": 10000,
        "read1_adapter_sequence": "AGATCGGAAGAGC",
        "read2_adapter_sequence": "AGATCGGAAGAGC",
        "read1_adapter_counts": {"AGATCGGAAGAGC": 300},
        "read2_adapter_counts": {"AGATCGGAAGAGC": 280}
    });
    let result = merge_adapter_cutting(&[a1, a2]);
    assert_eq!(result["adapter_trimmed_reads"].as_i64(), Some(600));
    assert_eq!(result["adapter_trimmed_bases"].as_i64(), Some(15_000));
    assert_eq!(result["read1_adapter_sequence"].as_str(), Some("AGATCGGAAGAGC"));
    assert_eq!(
        result["read1_adapter_counts"]["AGATCGGAAGAGC"].as_i64(),
        Some(450)
    );
    assert_eq!(
        result["read2_adapter_counts"]["AGATCGGAAGAGC"].as_i64(),
        Some(420)
    );
}

#[test]
fn test_merge_empty_list() {
    let result = merge_fastp_jsons(&[]);
    assert!(result.as_object().unwrap().is_empty());
}

// ── Integration tests using fixture files ─────────────────────────────────────

#[test]
fn test_merge_two_fixture_files_summary_totals() {
    let c1 = load_fixture("fastp_chunk1.json");
    let c2 = load_fixture("fastp_chunk2.json");
    let merged = merge_fastp_jsons(&[c1, c2]);

    let before = &merged["summary"]["before_filtering"];
    assert_eq!(before["total_reads"].as_i64(), Some(3000));
    assert_eq!(before["total_bases"].as_i64(), Some(450_000));
    assert_eq!(before["q20_bases"].as_i64(), Some(420_000));
    assert_eq!(before["q30_bases"].as_i64(), Some(360_000));
}

#[test]
fn test_merge_two_fixture_files_rates_recalculated() {
    let c1 = load_fixture("fastp_chunk1.json");
    let c2 = load_fixture("fastp_chunk2.json");
    let merged = merge_fastp_jsons(&[c1, c2]);

    let before = &merged["summary"]["before_filtering"];
    assert_approx_eq(
        before["q20_rate"].as_f64().unwrap(),
        420_000.0 / 450_000.0,
        "q20_rate",
    );
    assert_approx_eq(
        before["q30_rate"].as_f64().unwrap(),
        360_000.0 / 450_000.0,
        "q30_rate",
    );
    // gc_content: (0.42*150000 + 0.40*300000) / 450000 = 183000/450000
    assert_approx_eq(
        before["gc_content"].as_f64().unwrap(),
        183_000.0 / 450_000.0,
        "gc_content",
    );
    assert_eq!(before["read1_mean_length"].as_i64(), Some(150));
}

#[test]
fn test_merge_two_fixture_files_duplication_weighted() {
    let c1 = load_fixture("fastp_chunk1.json");
    let c2 = load_fixture("fastp_chunk2.json");
    let merged = merge_fastp_jsons(&[c1, c2]);

    // rate = (0.1*1000 + 0.2*2000) / 3000 = 500/3000
    assert_approx_eq(
        merged["duplication"]["rate"].as_f64().unwrap(),
        500.0 / 3000.0,
        "duplication rate",
    );
}

#[test]
fn test_merge_two_fixture_files_insert_size() {
    let c1 = load_fixture("fastp_chunk1.json");
    let c2 = load_fixture("fastp_chunk2.json");
    let merged = merge_fastp_jsons(&[c1, c2]);

    let insert = &merged["insert_size"];
    assert_eq!(insert["unknown"].as_i64(), Some(150));
    assert_eq!(insert["peak"].as_u64(), Some(7));
    let hist = insert["histogram"].as_array().unwrap();
    assert_eq!(hist[1].as_i64(), Some(15));
    assert_eq!(hist[7].as_i64(), Some(600));
    assert_eq!(hist[10].as_i64(), Some(150));
}

#[test]
fn test_merge_two_fixture_files_filtering_result() {
    let c1 = load_fixture("fastp_chunk1.json");
    let c2 = load_fixture("fastp_chunk2.json");
    let merged = merge_fastp_jsons(&[c1, c2]);

    let fr = &merged["filtering_result"];
    assert_eq!(fr["passed_filter_reads"].as_i64(), Some(2850));
    assert_eq!(fr["low_quality_reads"].as_i64(), Some(90));
    assert_eq!(fr["too_many_N_reads"].as_i64(), Some(30));
    assert_eq!(fr["too_short_reads"].as_i64(), Some(30));
    assert_eq!(fr["too_long_reads"].as_i64(), Some(0));
}

#[test]
fn test_merge_two_fixture_files_adapter_cutting() {
    let c1 = load_fixture("fastp_chunk1.json");
    let c2 = load_fixture("fastp_chunk2.json");
    let merged = merge_fastp_jsons(&[c1, c2]);

    let ac = &merged["adapter_cutting"];
    assert_eq!(ac["adapter_trimmed_reads"].as_i64(), Some(600));
    assert_eq!(ac["adapter_trimmed_bases"].as_i64(), Some(15_000));
    assert_eq!(ac["read1_adapter_sequence"].as_str(), Some("AGATCGGAAGAGC"));
    assert_eq!(
        ac["read1_adapter_counts"]["AGATCGGAAGAGC"].as_i64(),
        Some(450)
    );
    assert_eq!(
        ac["read2_adapter_counts"]["AGATCGGAAGAGC"].as_i64(),
        Some(420)
    );
}

#[test]
fn test_merge_two_fixture_files_per_read_sections() {
    let c1 = load_fixture("fastp_chunk1.json");
    let c2 = load_fixture("fastp_chunk2.json");
    let merged = merge_fastp_jsons(&[c1, c2]);

    let r1_before = &merged["read1_before_filtering"];
    assert_eq!(r1_before["total_reads"].as_i64(), Some(1500));
    assert_eq!(r1_before["total_bases"].as_i64(), Some(225_000));
    assert_eq!(r1_before["q20_bases"].as_i64(), Some(210_000));
    // total_cycles preserved from first chunk
    assert_eq!(r1_before["total_cycles"].as_i64(), Some(150));
    // kmer counts summed
    assert_eq!(r1_before["kmer_count"]["AAAAA"].as_i64(), Some(300));
    assert_eq!(r1_before["kmer_count"]["AAAAC"].as_i64(), Some(150));
    // quality curves averaged (both chunks have same values, so average = same)
    let mean = r1_before["quality_curves"]["mean"].as_array().unwrap();
    assert_approx_eq(mean[0].as_f64().unwrap(), 35.0, "quality mean[0]");

    let r1_after = &merged["read1_after_filtering"];
    assert_eq!(r1_after["total_reads"].as_i64(), Some(1425));
    assert_eq!(r1_after["kmer_count"]["AAAAA"].as_i64(), Some(285));

    let r2_before = &merged["read2_before_filtering"];
    assert_eq!(r2_before["total_reads"].as_i64(), Some(1500));

    let r2_after = &merged["read2_after_filtering"];
    assert_eq!(r2_after["total_reads"].as_i64(), Some(1425));
}

#[test]
fn test_merge_preserves_fastp_version_from_first_file() {
    let c1 = load_fixture("fastp_chunk1.json");
    let c2 = load_fixture("fastp_chunk2.json");
    let merged = merge_fastp_jsons(&[c1, c2]);
    assert_eq!(
        merged["summary"]["fastp_version"].as_str(),
        Some("0.23.2")
    );
}

#[test]
fn test_merge_single_file_is_identity_for_counts() {
    let c1 = load_fixture("fastp_chunk1.json");
    let merged = merge_fastp_jsons(&[c1]);

    let before = &merged["summary"]["before_filtering"];
    assert_eq!(before["total_reads"].as_i64(), Some(1000));
    assert_eq!(before["total_bases"].as_i64(), Some(150_000));
    // q20_rate should be recalculated from totals (should match original)
    assert_approx_eq(
        before["q20_rate"].as_f64().unwrap(),
        140_000.0 / 150_000.0,
        "q20_rate single file",
    );
    assert_eq!(merged["duplication"]["rate"].as_f64(), Some(0.1));
}
