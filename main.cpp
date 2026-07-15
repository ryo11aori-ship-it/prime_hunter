#pragma GCC optimize("O3,unroll-loops")
#include <iostream>
#include <vector>
#include <chrono>
#include <gmp.h>

using namespace std;
using namespace std::chrono;

// q <= 3581（500個目の奇素数）なので、乗算は絶対に64bit整数でオーバーフローしない。
// 重い __int128 キャストを排除して極限まで軽量化
inline long long power_mod(long long base, long long exp, long long mod) {
    long long res = 1;
    base %= mod;
    if (base == 0) return 0;
    while (exp > 0) {
        if (exp & 1) res = (res * base) % mod;
        base = (base * base) % mod;
        exp >>= 1;
    }
    return res;
}

inline long long p_pow_p_mod_q(long long p, long long q) {
    if (p % q == 0) return 0;
    long long exp = p % (q - 1);
    return power_mod(p, exp, q);
}

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
    const int TIME_LIMIT_SEC = 290; // 5分未満（290秒）で確実に安全停止

    // 探索上限を「25,000」（素数2,762個）に調整。最大桁数は約11万桁。
    // これならすり抜けたGMP計算も1回あたり0.1秒未満で終わるため、ハングしません。
    int limit = 25000; 
    cout << "Generating primes up to " << limit << "..." << endl;
    vector<int> primes = generate_primes(limit);
    int num_primes = primes.size();
    cout << "Generated " << num_primes << " primes." << endl;

    // フィルター用奇素数を500個（3から3581まで）使用
    // メルテンスの定理より、すり抜け率は約 0.56 / ln(3581) ≒ 6.8%
    vector<int> filter_primes;
    for (int p : primes) {
        if (p == 2) continue;
        filter_primes.push_back(p);
        if (filter_primes.size() >= 500) break;
    }
    int F = filter_primes.size();
    cout << "Using " << F << " filter primes (from " << filter_primes.front() << " to " << filter_primes.back() << ")." << endl;

    // ローリングバッファ
    vector<vector<long long>> rem(3, vector<long long>(F));
    for (int j = 0; j < F; j++) {
        rem[0][j] = p_pow_p_mod_q(primes[0], filter_primes[j]);
        rem[1][j] = p_pow_p_mod_q(primes[1], filter_primes[j]);
    }

    int found_count = 0;
    long long checked_count = 0;
    cout << "Starting safe and fast search..." << endl;

    for (int i = 2; i < num_primes; i++) {
        // 1000ループごとに時間チェック
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

        for (int j = 0; j < F; j++) {
            rem[c_idx][j] = p_pow_p_mod_q(c, filter_primes[j]);
        }

        bool pass_filter = true;
        for (int j = 0; j < F; j++) {
            long long q = filter_primes[j];
            long long sum_rem = (rem[a_idx][j] + rem[b_idx][j] + rem[c_idx][j]) % q;

            if (sum_rem == 0) {
                pass_filter = false;
                break;
            }
        }

        if (pass_filter) {
            // GMP計算に突入するまさにその直前でタイムアウトを厳重チェック
            auto elapsed = duration_cast<seconds>(steady_clock::now() - start_time).count();
            if (elapsed >= TIME_LIMIT_SEC) {
                cout << "\n[Timeout] " << TIME_LIMIT_SEC << " seconds reached inside hot path. Stopping." << endl;
                break;
            }

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
