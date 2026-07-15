#pragma GCC optimize("O3,unroll-loops")
#include <iostream>
#include <vector>
#include <chrono>
#include <gmp.h>

using namespace std;
using namespace std::chrono;

// 64ビット高速べき乗剰余 (base^exp % mod)
long long power_mod(long long base, long long exp, long long mod) {
    long long res = 1;
    base %= mod;
    if (base == 0) return 0;
    while (exp > 0) {
        if (exp & 1) res = (__int128)res * base % mod; // オーバーフロー防止用の128bitキャスト
        base = (__int128)base * base % mod;
        exp >>= 1;
    }
    return res;
}

// p^p % q を高速計算
long long p_pow_p_mod_q(long long p, long long q) {
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
    const int TIME_LIMIT_SEC = 300; // 5分

    // 探索上限を「5000万（50M）」まで超絶拡張（素数は約300万個存在）
    int limit = 50000000; 
    cout << "Generating primes up to " << limit << "..." << endl;
    vector<int> primes = generate_primes(limit);
    int num_primes = primes.size();
    cout << "Generated " << num_primes << " primes." << endl;

    // フィルター用素数として、最初の200個の素数（3から1223まで）を使用
    int num_filter_primes = min(200, num_primes);
    vector<int> filter_primes(primes.begin(), primes.begin() + num_filter_primes);

    int found_count = 0;
    long long checked_count = 0;
    cout << "Starting hyper-optimized search with 200-stage prime filter..." << endl;

    for (int i = 0; i < num_primes - 2; i++) {
        // 10000ループごとにタイムアウトを判定（オーバーヘッドを極限まで削減）
        if (i % 10000 == 0) {
            auto elapsed = duration_cast<seconds>(steady_clock::now() - start_time).count();
            if (elapsed >= TIME_LIMIT_SEC) {
                cout << "\n[Timeout] 5 minutes reached. Stopping search." << endl;
                break;
            }
        }

        int a = primes[i];
        int b = primes[i + 1];
        int c = primes[i + 2];
        checked_count++;

        // 200段階の高速Moduloフィルター
        bool pass_filter = true;
        for (int q : filter_primes) {
            long long rem = 0;
            rem = (rem + p_pow_p_mod_q(a, q)) % q;
            rem = (rem + p_pow_p_mod_q(b, q)) % q;
            rem = (rem + p_pow_p_mod_q(c, q)) % q;

            if (rem == 0) {
                pass_filter = false;
                break; // 合成数確定なので即座に次の3つ組へ
            }
        }

        // 奇跡的に200段階のフィルターをすべて突破した「本物の候補」のみ、GMPで実計算
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
