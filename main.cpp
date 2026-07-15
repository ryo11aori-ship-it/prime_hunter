#pragma GCC optimize("O3,unroll-loops")
#include <iostream>
#include <vector>
#include <chrono>
#include <gmp.h>

using namespace std;
using namespace std::chrono;

// 64ビット高速べき乗剰余 (base^exp % mod)
inline long long power_mod(long long base, long long exp, long long mod) {
    long long res = 1;
    base %= mod;
    if (base == 0) return 0;
    while (exp > 0) {
        if (exp & 1) res = (long long)((__int128)res * base % mod); // 128bitでオーバーフロー防止
        base = (long long)((__int128)base * base % mod);
        exp >>= 1;
    }
    return res;
}

// p^p % q を高速計算
inline long long p_pow_p_mod_q(long long p, long long q) {
    if (p % q == 0) return 0;
    long long exp = p % (q - 1);
    return power_mod(p, exp, q);
}

// エラトステネスの篩による高速素数生成
vector<int> generate_primes(int limit) {
    vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int p = 2; p * p <= limit; p++) {
        if (is_prime[p]) {
            for (int i = p * p; i <= limit; i += p)
                is_prime[i] = false;
        }
    }
    vector<int> primes;
    for (int p = 2; p <= limit; p++) {
        if (is_prime[p]) primes.push_back(p);
    }
    return primes;
}

int main() {
    setbuf(stdout, NULL);
    auto start_time = steady_clock::now();
    const int TIME_LIMIT_SEC = 290; // 安全のため290秒で終了

    // 探索上限を 50万 に設定（素数は38,000個以上。これでも従来の数万倍広い）
    int limit = 500000; 
    cout << "Generating primes up to " << limit << "..." << endl;
    vector<int> primes = generate_primes(limit);
    int num_primes = primes.size();
    cout << "Generated " << num_primes << " primes." << endl;

    // フィルター用奇素数を1000個用意（3から7927までの素数）
    vector<int> filter_primes;
    for (int p : primes) {
        if (p == 2) continue;
        filter_primes.push_back(p);
        if (filter_primes.size() >= 1000) break;
    }
    int F = filter_primes.size();
    cout << "Using " << F << " filter primes (from " << filter_primes.front() << " to " << filter_primes.back() << ")." << endl;

    // ローリングバッファ用の剰余テーブル [3][F]
    // rem[i % 3][j] に primes[i]^primes[i] % filter_primes[j] を格納
    vector<vector<long long>> rem(3, vector<long long>(F));

    // 最初の2つの素数(2, 3)の剰余をあらかじめ計算しておく
    for (int j = 0; j < F; j++) {
        rem[0][j] = p_pow_p_mod_q(primes[0], filter_primes[j]);
        rem[1][j] = p_pow_p_mod_q(primes[1], filter_primes[j]);
    }

    int found_count = 0;
    long long checked_count = 0;
    cout << "Starting hyper-optimized search..." << endl;

    for (int i = 2; i < num_primes; i++) {
        // 1000ループごとにタイムアウト判定（ラグをミリ秒以下に）
        if (i % 1000 == 0) {
            auto elapsed = duration_cast<seconds>(steady_clock::now() - start_time).count();
            if (elapsed >= TIME_LIMIT_SEC) {
                cout << "\n[Timeout] " << TIME_LIMIT_SEC << " seconds reached. Stopping search." << endl;
                break;
            }
        }

        int a = primes[i - 2];
        int b = primes[i - 1];
        int c = primes[i];
        checked_count++;

        int c_idx = i % 3;
        int a_idx = (i - 2) % 3;
        int b_idx = (i - 1) % 3;

        // 新しく入ってきた c^c % q のみを計算（重複計算の排除）
        for (int j = 0; j < F; j++) {
            rem[c_idx][j] = p_pow_p_mod_q(c, filter_primes[j]);
        }

        // 1000段階の高速Moduloフィルター
        bool pass_filter = true;
        for (int j = 0; j < F; j++) {
            long long q = filter_primes[j];
            long long sum_rem = (rem[a_idx][j] + rem[b_idx][j] + rem[c_idx][j]) % q;

            if (sum_rem == 0) {
                pass_filter = false;
                break; // 合成数確定なので即座にスキップ
            }
        }

        // 奇跡的に1000個のフィルターをすべて突破した「本物の候補」のみ、GMPで実計算
        if (pass_filter) {
            mpz_t gmp_a, gmp_b, gmp_c, gmp_sum;
            mpz_inits(gmp_a, gmp_b, gmp_c, gmp_sum, NULL);

            mpz_ui_pow_ui(gmp_a, a, a);
            mpz_ui_pow_ui(gmp_b, b, b);
            mpz_ui_pow_ui(gmp_c, c, c);

            mpz_add(gmp_sum, gmp_a, gmp_b);
            mpz_add(gmp_sum, gmp_sum, gmp_c);

            // ミラー・ラビン素数判定
            if (mpz_probab_prime_p(gmp_sum, 25) > 0) {
                found_count++;
                int digits = mpz_sizeinbase(gmp_sum, 10);
                cout << "Found! a=" << a << ", b=" << b << ", c=" << c 
                     << " (Sum has " << digits << " digits)" << endl;
            }

            mpz_clears(gmp_a, gmp_b, gmp_c, gmp_sum, NULL);
        }
    }

    auto end_time = steady_clock::now();
    auto total_elapsed = duration_cast<milliseconds>(end_time - start_time).count() / 1000.0;
    cout << "Search finished." << endl;
    cout << "Total checked triples: " << checked_count << endl;
    cout << "Total found: " << found_count << endl;
    cout << "Execution time: " << total_elapsed << " seconds" << endl;

    return 0;
}
