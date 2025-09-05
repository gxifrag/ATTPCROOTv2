void merge_selected_runs() {
    const char* root_files[] = {
        "run_0024_2H.root", "run_0022_2H.root", "run_0115_2H.root", "run_0066_2H.root",
        "run_0124_2H.root", "run_0102_2H.root", "run_0103_2H.root", "run_0044_2H.root",
        "run_0051_2H.root", "run_0110_2H.root", "run_0073_2H.root", "run_0123_2H.root",
        "run_0016_2H.root", "run_0089_2H.root", "run_0017_2H.root", "run_0062_2H.root",
        "run_0114_2H.root", "run_0128_2H.root", "run_0014_2H.root", "run_0119_2H.root",
        "run_0050_2H.root", "run_0069_2H.root", "run_0096_2H.root", "run_0105_2H.root",
        "run_0068_2H.root", "run_0107_2H.root", "run_0070_2H.root", "run_0034_2H.root",
        "run_0023_2H.root", "run_0130_2H.root", "run_0036_2H.root", "run_0101_2H.root",
        "run_0106_2H.root", "run_0060_2H.root", "run_0071_2H.root", "run_0061_2H.root",
        "run_0084_2H.root", "run_0113_2H.root", "run_0075_2H.root", "run_0043_2H.root",
        "run_0021_2H.root", "run_0133_2H.root", "run_0104_2H.root", "run_0129_2H.root",
        "run_0088_2H.root", "run_0042_2H.root", "run_0117_2H.root", "run_0037_2H.root",
        "run_0013_2H.root", "run_0030_2H.root", "run_0100_2H.root", "run_0015_2H.root",
        "run_0122_2H.root", "run_0099_2H.root", "run_0074_2H.root", "run_0033_2H.root",
        "run_0020_2H.root", "run_0028_2H.root", "run_0097_2H.root", "run_0029_2H.root",
        "run_0121_2H.root", "run_0120_2H.root", "run_0127_2H.root", "run_0072_2H.root",
        "run_0083_2H.root", "run_0056_2H.root", "run_0131_2H.root", "run_0085_2H.root",
        "run_0054_2H.root"
    };

    const int n_files = sizeof(root_files) / sizeof(root_files[0]);
    TString treeName = "parquettree"; // Change if your TTree is named differently
    TChain chain(treeName);

    // Add full path to each file
    TString base_path = "/home/georgina/ATTPC/rootfiles/";
    for (int i = 0; i < n_files; ++i) {
        std::cout << "Adding file: " << root_files[i] << std::endl;
        chain.Add(base_path + root_files[i]);
    }

    TFile* outfile = TFile::Open("merged_results.root", "RECREATE");
    chain.Merge(outfile, 0, "keep"); // keep = preserve branches as-is
    outfile->Close();

    std::cout << "✅ Merge complete. Output saved to merged_results.root" << std::endl;
}
