#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <zlib.h>
#include <omp.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_cdf.h>

using namespace std;

struct CSVData {
    vector<string> colNames;
    vector<string> rowNames;
    vector<vector<double>> data;
};

// Write string to gzipped file
void writeToGzip(gzFile file, const string& str) {
    gzwrite(file, str.c_str(), str.length());
}

// Parse CSV file (write.csv2 format: semicolon separator, comma decimal)
// Handles both regular and gzipped files
CSVData parseCSV(const string& filename) {
    CSVData csv;
    bool isGzipped = false;

    // Check if file is gzipped by extension
    if (filename.size() > 3 && filename.substr(filename.size() - 3) == ".gz") {
        isGzipped = true;
    }

    vector<string> lines;

    if (isGzipped) {
        // Read gzipped file
        gzFile gzfile = gzopen(filename.c_str(), "rb");
        if (!gzfile) {
            cerr << "Error: Cannot open gzipped file " << filename << endl;
            exit(1);
        }

        char buffer[65536];
        string currentLine = "";

        while (true) {
            int bytesRead = gzread(gzfile, buffer, sizeof(buffer) - 1);
            if (bytesRead < 0) {
                cerr << "Error reading gzipped file " << filename << endl;
                gzclose(gzfile);
                exit(1);
            }
            if (bytesRead == 0) break;

            buffer[bytesRead] = '\0';

            for (int i = 0; i < bytesRead; i++) {
                if (buffer[i] == '\n') {
                    lines.push_back(currentLine);
                    currentLine = "";
                } else if (buffer[i] != '\r') {
                    currentLine += buffer[i];
                }
            }
        }

        if (!currentLine.empty()) {
            lines.push_back(currentLine);
        }

        gzclose(gzfile);
    } else {
        // Read regular file
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error: Cannot open file " << filename << endl;
            exit(1);
        }

        string line;
        while (getline(file, line)) {
            lines.push_back(line);
        }
        file.close();
    }

    // Process lines
    bool firstRow = true;
    for (const string& line : lines) {
        stringstream ss(line);
        string cell;
        vector<string> row;

        while (getline(ss, cell, ';')) {
            // Remove quotes if present
            if (!cell.empty() && cell.front() == '"' && cell.back() == '"') {
                cell = cell.substr(1, cell.length() - 2);
            }
            row.push_back(cell);
        }
        if (firstRow) {
            // First row contains column names (skip first cell which is empty/"")
            for (size_t i = 1; i < row.size(); i++) {
                csv.colNames.push_back(row[i]);
            }
            firstRow = false;
        } else {
            // First cell is row name, rest are data
            if (!row.empty()) {
                csv.rowNames.push_back(row[0]);
                vector<double> dataRow;
                for (size_t i = 1; i < row.size(); i++) {
                    // Replace comma with dot for decimal point
                    string val = row[i];
                    replace(val.begin(), val.end(), ',', '.');
                    try {
                        double value = stod(val);
                        // Check for NaN or Inf
                        if (!isfinite(value)) {
                            cerr << "Warning: Non-finite value found, replacing with 0" << endl;
                            value = 0.0;
                        }
                        dataRow.push_back(value);
                    } catch (...) {
                        cerr << "Error parsing value: " << row[i] << ", using 0.0" << endl;
                        dataRow.push_back(0.0);
                    }
                }
                if (!dataRow.empty()) {
                    csv.data.push_back(dataRow);
                }
            }
        }
    }

    if (csv.data.empty() || csv.colNames.empty()) {
        cerr << "Error: No valid data found in CSV file" << endl;
        exit(1);
    }

    return csv;
}
// Incomplete beta function (needed for t-distribution)
double incompleteBeta(double x, double a, double b) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;

    // Use continued fraction approximation
    double bt = exp(lgamma(a + b) - lgamma(a) - lgamma(b) +
                    a * log(x) + b * log(1.0 - x));

    if (x < (a + 1.0) / (a + b + 2.0)) {
        // Use continued fraction directly
        double qab = a + b;
        double qap = a + 1.0;
        double qam = a - 1.0;
        double c = 1.0;
        double d = 1.0 - qab * x / qap;

        if (fabs(d) < 1e-30) d = 1e-30;
        d = 1.0 / d;
        double h = d;

        for (int m = 1; m <= 100; m++) {
            int m2 = 2 * m;
            double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1.0 + aa * d;
            if (fabs(d) < 1e-30) d = 1e-30;
            c = 1.0 + aa / c;
            if (fabs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d;
            h *= d * c;

            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = 1.0 + aa * d;
            if (fabs(d) < 1e-30) d = 1e-30;
            c = 1.0 + aa / c;
            if (fabs(c) < 1e-30) c = 1e-30;
            d = 1.0 / d;
            double del = d * c;
            h *= del;

            if (fabs(del - 1.0) < 1e-10) break;
        }

        return bt * h / a;
    } else {
        // Use symmetry relation
        return 1.0 - incompleteBeta(1.0 - x, b, a);
    }
}
// Calculate p-value for correlation coefficient using t-distribution
double correlationPValue(double r, int n) {
    if (n < 3) return 1.0;

    // Handle perfect correlations
    if (fabs(r) >= 0.9999999) {
        return (fabs(r) >= 1.0) ? 0.0 : 1e-300;
    }

    // Calculate t-statistic
    double t = r * sqrt(n - 2.0) / sqrt(1.0 - r * r);
    int df = n - 2;

    // Calculate two-tailed p-value using incomplete beta function
    // P(|T| > |t|) where T ~ t(df)
    double x = df / (df + t * t);
    double pval = incompleteBeta(x, df / 2.0, 0.5);

    // Ensure p-value is in valid range
    if (pval < 0.0) pval = 0.0;
    if (pval > 1.0) pval = 1.0;

    return pval;
}

// Get column from data matrix
vector<double> getColumn(const vector<vector<double>>& data, int col) {
    vector<double> column;
    column.reserve(data.size());
    for (const auto& row : data) {
        if (col < (int)row.size()) {
            column.push_back(row[col]);
        }
    }
    return column;
}

// Calculate Pearson correlation using GSL
double pearsonCorrelation(const vector<double>& x, const vector<double>& y) {
    if (x.size() != y.size() || x.empty()) return 0.0;
    if (x.size() < 2) return 0.0;

    // Check for constant vectors
    bool x_constant = true, y_constant = true;
    for (size_t i = 1; i < x.size(); i++) {
        if (x[i] != x[0]) x_constant = false;
        if (y[i] != y[0]) y_constant = false;
    }
    if (x_constant || y_constant) return 0.0;

    return gsl_stats_correlation(x.data(), 1, y.data(), 1, x.size());
}
// Calculate Spearman correlation using GSL
double spearmanCorrelation(const vector<double>& x, const vector<double>& y) {
    if (x.size() != y.size() || x.empty()) return 0.0;
    if (x.size() < 2) return 0.0;

    size_t n = x.size();

    // Check for constant vectors
    bool x_constant = true, y_constant = true;
    for (size_t i = 1; i < n; i++) {
        if (x[i] != x[0]) x_constant = false;
        if (y[i] != y[0]) y_constant = false;
    }
    if (x_constant || y_constant) return 0.0;

    // Create copies to avoid modifying input
    vector<double> x_copy = x;
    vector<double> y_copy = y;

    // Allocate work space (must be 2*n for GSL)
    vector<double> work(2 * n);

    double result = gsl_stats_spearman(x_copy.data(), 1, y_copy.data(), 1, n, work.data());

    return result;
}

// Mode 1: All pairwise correlations within one CSV
void mode1(const string& filename, const string& outputFile = "", int numThreads = 1) {
    CSVData csv = parseCSV(filename);
    int nCols = csv.colNames.size();

    // Calculate total number of pairs
    long long totalPairs = ((long long)nCols * (nCols - 1)) / 2;

    // Memory usage estimate
    size_t estimatedMemoryMB = (csv.data.size() * nCols * sizeof(double)) / (1024 * 1024);

    if (outputFile.empty()) {
        cout << "\n=== Mode 1: Within-matrix correlations ===" << endl;
        cout << "File: " << filename << endl;
        cout << "Columns: " << nCols << ", Rows: " << csv.data.size() << endl;
        cout << "Total pairs: " << totalPairs << endl;
        cout << "Estimated memory usage: ~" << estimatedMemoryMB << " MB" << endl;
        cout << "Threads: " << numThreads << endl << endl;

        cout << left << setw(25) << "Column 1"
             << setw(25) << "Column 2"
             << setw(12) << "Pearson"
             << setw(12) << "P-value"
             << setw(12) << "Spearman"
             << setw(12) << "P-value" << endl;
        cout << string(98, '-') << endl;
    }
    // Set number of threads
    omp_set_num_threads(numThreads);

    // Use file output or vector storage based on output mode
    gzFile outfile = NULL;
    if (!outputFile.empty()) {
        outfile = gzopen(outputFile.c_str(), "wb");
        if (!outfile) {
            cerr << "Error: Cannot create output file " << outputFile << endl;
            exit(1);
        }
        writeToGzip(outfile, "Column1\tColumn2\tPearson\tPearson_pval\tSpearman\tSpearman_pval\n");
    }

    // Process in parallel, but write sequentially to avoid corruption
    #pragma omp parallel for schedule(dynamic) ordered
    for (long long idx = 0; idx < totalPairs; idx++) {
        // Convert linear index to (i, j) pair in O(1) via closed-form triangular inversion.
        // Total pairs before row i = i*(2*nCols - i - 1)/2, so solving for i gives:
        int i = (int)((2.0*nCols - 1.0 - sqrt((2.0*nCols - 1.0)*(2.0*nCols - 1.0) - 8.0*idx)) / 2.0);
        // Correct for potential floating-point rounding
        if ((long long)(i + 1) * (2 * nCols - i - 2) / 2 <= idx) i++;
        long long pairs_before_i = (long long)i * (2 * nCols - i - 1) / 2;
        int j = (int)(idx - pairs_before_i + i + 1);

        vector<double> col1 = getColumn(csv.data, i);
        vector<double> col2 = getColumn(csv.data, j);

        int n = col1.size();
        double pearson = pearsonCorrelation(col1, col2);
        double spearman = spearmanCorrelation(col1, col2);
        double pearsonPval = correlationPValue(pearson, n);
        double spearmanPval = correlationPValue(spearman, n);
        #pragma omp ordered
        {
            if (outputFile.empty()) {
                cout << left << setw(25) << csv.colNames[i]
                     << setw(25) << csv.colNames[j]
                     << setw(12) << fixed << setprecision(6) << pearson
                     << setw(12) << scientific << setprecision(3) << pearsonPval
                     << setw(12) << fixed << setprecision(6) << spearman
                     << setw(12) << scientific << setprecision(3) << spearmanPval << endl;
            } else {
                stringstream ss;
                ss << csv.colNames[i] << "\t" << csv.colNames[j] << "\t"
                   << fixed << setprecision(6) << pearson << "\t"
                   << scientific << setprecision(6) << pearsonPval << "\t"
                   << fixed << setprecision(6) << spearman << "\t"
                   << scientific << setprecision(6) << spearmanPval << "\n";
                writeToGzip(outfile, ss.str());
            }
        }
    }

    if (outfile) {
        gzclose(outfile);
        cout << "Results saved to " << outputFile << endl;
    }
}

// Mode 2: Cross-matrix correlations
void mode2(const string& file1, const string& file2, const string& outputFile = "", int numThreads = 1) {
    CSVData csv1 = parseCSV(file1);
    CSVData csv2 = parseCSV(file2);

    if (csv1.colNames.empty() || csv2.colNames.empty()) {
        cerr << "Error: One or both files have no columns" << endl;
        exit(1);
    }

    if (csv1.data.size() != csv2.data.size()) {
        cerr << "Warning: Files have different number of rows. Using minimum." << endl;
    }

    size_t minRows = min(csv1.data.size(), csv2.data.size());

    if (minRows < 3) {
        cerr << "Error: Need at least 3 rows for correlation analysis" << endl;
        exit(1);
    }

    long long totalPairs = (long long)csv1.colNames.size() * csv2.colNames.size();
    // Memory usage estimate
    size_t estimatedMemoryMB = ((csv1.data.size() * csv1.colNames.size() +
                                  csv2.data.size() * csv2.colNames.size()) * sizeof(double)) / (1024 * 1024);

    if (outputFile.empty()) {
        cout << "\n=== Mode 2: Cross-matrix correlations ===" << endl;
        cout << "File 1: " << file1 << " (" << csv1.colNames.size() << " columns, "
             << csv1.data.size() << " rows)" << endl;
        cout << "File 2: " << file2 << " (" << csv2.colNames.size() << " columns, "
             << csv2.data.size() << " rows)" << endl;
        cout << "Using " << minRows << " rows for correlation" << endl;
        cout << "Total pairs: " << totalPairs << endl;
        cout << "Estimated memory usage: ~" << estimatedMemoryMB << " MB" << endl;
        cout << "Threads: " << numThreads << endl << endl;

        cout << left << setw(25) << "Column (File 1)"
             << setw(25) << "Column (File 2)"
             << setw(12) << "Pearson"
             << setw(12) << "P-value"
             << setw(12) << "Spearman"
             << setw(12) << "P-value" << endl;
        cout << string(98, '-') << endl;
    }

    // Set number of threads
    omp_set_num_threads(numThreads);

    // Use file output or direct printing to avoid storing all results
    gzFile outfile = NULL;
    if (!outputFile.empty()) {
        outfile = gzopen(outputFile.c_str(), "wb");
        if (!outfile) {
            cerr << "Error: Cannot create output file " << outputFile << endl;
            exit(1);
        }
        writeToGzip(outfile, "Column1\tColumn2\tPearson\tPearson_pval\tSpearman\tSpearman_pval\n");
    }

    // Parallelize the computation
    #pragma omp parallel for schedule(dynamic) ordered
    for (long long idx = 0; idx < totalPairs; idx++) {
        int i = idx / csv2.colNames.size();
        int j = idx % csv2.colNames.size();

        vector<double> col1 = getColumn(csv1.data, i);
        vector<double> col2 = getColumn(csv2.data, j);

        // Trim to minimum size
        col1.resize(minRows);
        col2.resize(minRows);
        // Validate columns
        if (col1.size() != col2.size() || col1.size() < 3) {
            #pragma omp ordered
            {
                if (outputFile.empty()) {
                    cout << left << setw(25) << csv1.colNames[i]
                         << setw(25) << csv2.colNames[j]
                         << setw(12) << "NA"
                         << setw(12) << "NA"
                         << setw(12) << "NA"
                         << setw(12) << "NA" << endl;
                } else {
                    stringstream ss;
                    ss << csv1.colNames[i] << "\t" << csv2.colNames[j] << "\t"
                       << "NA\tNA\tNA\tNA\n";
                    writeToGzip(outfile, ss.str());
                }
            }
            continue;
        }

        int n = minRows;
        double pearson = pearsonCorrelation(col1, col2);
        double spearman = spearmanCorrelation(col1, col2);
        double pearsonPval = correlationPValue(pearson, n);
        double spearmanPval = correlationPValue(spearman, n);

        #pragma omp ordered
        {
            if (outputFile.empty()) {
                cout << left << setw(25) << csv1.colNames[i]
                     << setw(25) << csv2.colNames[j]
                     << setw(12) << fixed << setprecision(6) << pearson
                     << setw(12) << scientific << setprecision(3) << pearsonPval
                     << setw(12) << fixed << setprecision(6) << spearman
                     << setw(12) << scientific << setprecision(3) << spearmanPval << endl;
            } else {
                stringstream ss;
                ss << csv1.colNames[i] << "\t" << csv2.colNames[j] << "\t"
                   << fixed << setprecision(6) << pearson << "\t"
                   << scientific << setprecision(6) << pearsonPval << "\t"
                   << fixed << setprecision(6) << spearman << "\t"
                   << scientific << setprecision(6) << spearmanPval << "\n";
                writeToGzip(outfile, ss.str());
            }
        }
    }

    if (outfile) {
        gzclose(outfile);
        cout << "Results saved to " << outputFile << endl;
    }
}
int main(int argc, char* argv[]) {
    string outputFile = "";
    int numThreads = 1;
    vector<string> inputFiles;

    // Parse command line arguments
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            cout << "CSV Correlation Calculator" << endl;
            cout << "=========================" << endl << endl;
            cout << "Calculate Pearson and Spearman correlation coefficients for CSV data." << endl;
            cout << "Uses GNU Scientific Library (GSL) for accurate statistical calculations." << endl << endl;
            cout << "USAGE:" << endl;
            cout << "  Mode 1 (within-matrix correlations):" << endl;
            cout << "    " << argv[0] << " <csv_file> [options]" << endl << endl;
            cout << "  Mode 2 (cross-matrix correlations):" << endl;
            cout << "    " << argv[0] << " <csv_file1> <csv_file2> [options]" << endl << endl;
            cout << "OPTIONS:" << endl;
            cout << "  -o <file>      Save results to gzipped TSV file instead of stdout" << endl;
            cout << "  -t <N>         Use N threads for parallel computation (default: 1)" << endl;
            cout << "  --threads <N>  Same as -t" << endl;
            cout << "  -h, --help     Show this help message" << endl << endl;
            cout << "INPUT FILES:" << endl;
            cout << "  - CSV files with semicolon separator (R write.csv2 format)" << endl;
            cout << "  - First row contains column names, first column contains row names" << endl;
            cout << "  - Can be regular (.csv) or gzipped (.csv.gz) files" << endl;
            cout << "  - Decimal separator: comma (e.g., 3,14 for 3.14)" << endl << endl;
            cout << "OUTPUT:" << endl;
            cout << "  - Tab-separated values with columns:" << endl;
            cout << "    Column1, Column2, Pearson, Pearson_pval, Spearman, Spearman_pval" << endl;
            cout << "  - P-values are calculated using t-distribution (two-tailed test)" << endl;
            cout << "  - Printed to stdout or saved to gzipped file with -o option" << endl << endl;
            cout << "EXAMPLES:" << endl;
            cout << "  # Calculate within-matrix correlations" << endl;
            cout << "  " << argv[0] << " data.csv" << endl << endl;
            cout << "  # Use 8 threads and save to file" << endl;
            cout << "  " << argv[0] << " data.csv.gz -t 8 -o results.tsv.gz" << endl << endl;
            cout << "  # Calculate cross-matrix correlations with 4 threads" << endl;
            cout << "  " << argv[0] << " data1.csv data2.csv -t 4" << endl << endl;
            cout << "COMPILATION:" << endl;
            cout << "  g++ -o correlate correlate.cpp -std=c++11 -lz -fopenmp -lgsl -lgslcblas" << endl << endl;
            return 0;
        } else if (arg == "-o") {
            if (i + 1 < argc) {
                outputFile = argv[++i];
            } else {
                cerr << "Error: -o flag requires output filename" << endl;
                cerr << "Use -h or --help for usage information" << endl;
                return 1;
            }
        } else if (arg == "-t" || arg == "--threads") {
            if (i + 1 < argc) {
                numThreads = atoi(argv[++i]);
                if (numThreads < 1) {
                    cerr << "Error: Number of threads must be at least 1" << endl;
                    cerr << "Use -h or --help for usage information" << endl;
                    return 1;
                }
            } else {
                cerr << "Error: -t/--threads flag requires number of threads" << endl;
                cerr << "Use -h or --help for usage information" << endl;
                return 1;
            }
        } else {
            inputFiles.push_back(arg);
        }
    }

    if (inputFiles.empty() || inputFiles.size() > 2) {
        cerr << "Error: Invalid number of input files" << endl << endl;
        cout << "Usage:" << endl;
        cout << "  Mode 1: " << argv[0] << " <csv_file> [-o output.tsv.gz] [-t threads]" << endl;
        cout << "  Mode 2: " << argv[0] << " <csv_file1> <csv_file2> [-o output.tsv.gz] [-t threads]" << endl << endl;
        cout << "Use -h or --help for detailed usage information" << endl;
        return 1;
    }

    if (inputFiles.size() == 1) {
        mode1(inputFiles[0], outputFile, numThreads);
    } else {
        mode2(inputFiles[0], inputFiles[1], outputFile, numThreads);
    }

    return 0;
}
