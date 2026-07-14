#include<iostream>
#include<vector>
#include<gmp.h>
using namespace std;
bool is_prime(int n){
if(n<=1)return false;
if(n<=3)return true;
if(n%2==0||n%3==0)return false;
for(int i=5;i*i<=n;i+=6){
if(n%i==0||n%(i+2)==0)return false;
}
return true;
}
int main(){
vector<int> primes;
for(int i=3;i<2000;i+=2){
if(is_prime(i))primes.push_back(i);
}
int num_primes=primes.size();
vector<mpz_t> pows(num_primes);
for(int i=0;i<num_primes;i++){
mpz_init(pows[i]);
mpz_ui_pow_ui(pows[i],primes[i],primes[i]);
}
mpz_t sum;
mpz_init(sum);
int found_count=0;
for(int c=2;c<num_primes;c++){
for(int b=1;b<c;b++){
for(int a=0;a<b;a++){
mpz_add(sum,pows[a],pows[b]);
mpz_add(sum,sum,pows[c]);
if(mpz_probab_prime_p(sum,25)>0){
found_count++;
cout<<"["<<found_count<<"] a="<<primes[a]<<", b="<<primes[b]<<", c="<<primes[c]<<endl;
if(found_count>=100){
goto finish;
}
}
}
}
}
finish:
for(int i=0;i<num_primes;i++)mpz_clear(pows[i]);
mpz_clear(sum);
if(found_count<100)cout<<"Found only "<<found_count<<" in the given range."<<endl;
else cout<<"Successfully found 100 combinations."<<endl;
return 0;
}
