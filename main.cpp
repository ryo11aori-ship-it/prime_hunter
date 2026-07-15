#pragma GCC optimize("O3,unroll-loops")
#include <iostream>
#include <vector>
#include <chrono>
#include <gmp.h>

using namespace std;
using namespace std::chrono;

// エラトステネスの篩による高速な素数生成
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
    // バッファリングを無効化し、GitHub Actionsのログに即時反映させる
    setbuf(stdout, NULL);

    auto start_time = steady_clock::now();
    const int TIME_LIMIT_SEC = 300; // 5分でタイムアウト

    // 5分間のGMP累乗計算なら、素数は100万程度まで生成しておけば十分枯渇しない
    int limit = 1000000; 
    vector<int> primes = generate_primes(limit);
    int num_primes = primes.size();

    mpz_t p_pow[3];
    for (int i = 0; i < 3; i++) mpz_init(p_pow[i]);
    mpz_t sum;
    mpz_init(sum);

    // 最初の2つの素数(2, 3)の累乗だけ先に計算しておく
    mpz_ui_pow_ui(p_pow[0], primes[0], primes[0]);
    mpz_ui_pow_ui(p_pow[1], primes[1], primes[1]);

    int found_count = 0;
    cout << "Starting optimized search for 5 minutes..." << endl;

    for (int i = 0; i < num_primes - 2; i++) {
        // 100イテレーションごとにタイムアウトをチェック（clock呼び出しのオーバーヘッド削減）
        if (i % 100 == 0) {
            auto elapsed = duration_cast<seconds>(steady_clock::now() - start_time).count();
            if (elapsed >= TIME_LIMIT_SEC) {
                cout << "\n[Timeout] 5 minutes reached. Stopping search." << endl;
                break;
            }
        }

        int a = primes[i];
        int b = primes[i + 1];
        int c = primes[i + 2];

        // 剰余による数学的枝刈り (Modulo 3)
        // 奇素数pにおいて p^p ≡ p (mod 3)。ただしp=2のみ 2^2=4 ≡ 1 (mod 3)
        int mod_sum = 0;
        mod_sum += (a == 2) ? 1 : (a % 3);
        mod_sum += (b == 2) ? 1 : (b % 3);
        mod_sum += (c == 2) ? 1 : (c % 3);

        // a^a + b^b + c^c が3の倍数になる場合、素数ではない（自明な合成数）のでGMP計算をスキップ
        if (mod_sum % 3 != 0) {
            // 新しくウィンドウに入った c^c のみを計算し、古い a^a を上書き
            mpz_ui_pow_ui(p_pow[(i + 2) % 3], c, c);

            // 合計を計算
            mpz_add(sum, p_pow[0], p_pow[1]);
            mpz_add(sum, sum, p_pow[2]);

            // 確率的素数判定（ミラー・ラビン判定: 25回）
            if (mpz_probab_prime_p(sum, 25) > 0) {
                found_count++;
                cout << "Found! a=" << a << ", b=" << b << ", c=" << c << endl;
            }
        } else {
            // Modulo 3で弾かれた場合でも、ローリングバッファの整合性を保つために計算は必要
            mpz_ui_pow_ui(p_pow[(i + 2) % 3], c, c);
        }
    }

    cout << "Search finished. Total found: " << found_count << endl;

    // メモリ解放
    for (int i = 0; i < 3; i++) mpz_clear(p_pow[i]);
    mpz_clear(sum);

    return 0;
}
