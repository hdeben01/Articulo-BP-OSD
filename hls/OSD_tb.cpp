#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "OSD.h"


bool read_syndrome(std::istream& is, bit_t synd[N]) {
    int count = 0;
    std::string token;
    while (count < N && is >> token) {
        for (char c : token) {
            if (c == '0' || c == '1') {
                synd[count++] = (c - '0');
                if (count == N) break;
            }
        }
    }
    return (count == N);
}


bool read_probs(std::istream& is, value_t prob[M]) {
    int count = 0;
    std::string token;
    while (count < M && is >> token) {
        std::string clean = "";
        for (char c : token) {
            if (c != '[' && c != ']' && c != '*') {
                clean += c;
            }
        }
        if (!clean.empty()) {
            char* endptr = nullptr;
            double val = std::strtod(clean.c_str(), &endptr);
            if (endptr != clean.c_str()) {
                prob[count++] = (value_t)val;
            }
        }
    }
    return (count == M);
}


bool read_expected(std::istream& is, bit_t expected[M]) {
    int count = 0;
    std::string token;
    while (count < M && is >> token) {
        for (char c : token) {
            if (c == '0' || c == '1') {
                expected[count++] = (c - '0');
                if (count == M) break;
            }
        }
    }
    return (count == M);
}

int main() {
    std::ifstream file("input.txt");
    if (!file.is_open()) {
        std::cerr << "Error: Could not open input.txt" << std::endl;
        return 1;
    }

    bit_t synd[N];
    value_t prob[M];
    bit_t expected_sol[M];
    bit_t sol[M];

    int test_count = 0;
    int pass_count = 0;
    int total_bit_errors = 0;

    std::cout << "Starting OSD testbench with input.txt..." << std::endl;

    while (read_syndrome(file, synd)) {
        if (!read_probs(file, prob)) {
            std::cerr << "Error reading probabilities for testcase " << (test_count + 1) << std::endl;
            break;
        }
        if (!read_expected(file, expected_sol)) {
            std::cerr << "Error reading expected result for testcase " << (test_count + 1) << std::endl;
            break;
        }

        decode(synd, prob, sol);

        int bit_errors = 0;
        for (int j = 0; j < M; j++) {
            if (sol[j] != expected_sol[j]) {
                bit_errors++;
            }
        }

        test_count++;
        total_bit_errors += bit_errors;

        if (bit_errors == 0) {
            pass_count++;
            std::cout << "Test " << test_count << ": PASS" << std::endl;
        } else {
            std::cout << "Test " << test_count << ": FAIL (" << bit_errors << " bit errors)" << std::endl;
        }
    }

    file.close();

    std::cout << "\n========================================" << std::endl;
    std::cout << "Total tests executed: " << test_count << std::endl;
    std::cout << "Passed: " << pass_count << " / " << test_count << std::endl;
    std::cout << "Failed: " << (test_count - pass_count) << " / " << test_count << std::endl;
    std::cout << "Total bit errors: " << total_bit_errors << std::endl;
    std::cout << "========================================" << std::endl;

    if (pass_count == test_count && test_count > 0) {
        std::cout << "ALL TESTS PASSED!" << std::endl;
        return 0;
    } else {
        std::cout << "TESTBENCH FAILED!" << std::endl;
        return 1;
    }
}