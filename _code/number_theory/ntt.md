```cpp
const int mod = 16515073;
const int g = 5;

void reverse(int *d, int bits) {
	int n = (1 << bits), j = (1 << bits) >> 1;
	for (int i = 1; i < n - 1; i++) {
		if (i < j) swap(d[i], d[j]);
		int t = (1 << bits) >> 1;
		while (t <= j) {
			j -= t;
			t >>= 1;
		}
		j += t;
	}
}
void fft(int *d, int bits, int on) {
	reverse(d, bits);
	int n = (1 << bits);
	int aoe = g;
	if (on == -1) aoe = power(g, mod - 2);
	for (int s = 0; s < bits; s++) {
		int m = (1 << s) << 1, p = (1 << s);
		int w = power(aoe, (mod - 1) / m);
		for (int i = 0; i < n; i += m) {
			int cur = 1;
			for (int j = 0; j < p; j++) {
				int x = d[i + j], y = 1LL * d[i + p + j] * cur % mod;
				d[i + j] = (x + y) % mod;
				d[i + p + j] = (x - y) % mod;
				cur = 1LL * cur * w % mod;
			}
		}
	}
	if (on == -1) {
		int aoe = power(n, mod - 2);
		for (int i = 0; i < n; i++) d[i] = 1LL * d[i] * aoe % mod;
	}
}

int main() {
	int n, m;
	scanf("%d%d", &n, &m);
	for (int i = 0; i <= n; i++) scanf("%d", &a[i]);
	for (int i = 0; i <= m; i++) scanf("%d", &b[i]);
	
	int bits = 0;
	while ((1 << bits) < max(n + 1, m + 1)) bits++;
	bits++;
	
	int c = (1 << bits);
	fft(a, bits, 1);
	fft(b, bits, 1);
	for (int i = 0; i < c; i++) a[i] = 1LL * a[i] * b[i] % mod;
	fft(a, bits, -1);
	for (int i = 0; i <= n + m; i++) printf("%d ", (a[i] + mod) % mod);
	printf("\n");
	return 0;
}
```
