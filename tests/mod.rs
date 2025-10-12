// Removed obsolete test modules that used AssemblyGraphBuilder and other deprecated APIs:
// - comprehensive_test_suite
// - tdd_assembly_robustness_tests
// - test_edge_cases_bioinformatics
// - test_statistical_accuracy

pub mod test_assembly_pipeline;
pub mod test_core_data_structures;
pub mod test_graph_algorithms;
pub mod test_integration_end_to_end;
pub mod test_property_based_genomics;
pub mod test_synthetic_biological_data;
