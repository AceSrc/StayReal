```cpp
void exkmp(char *ch) {
	int n = strlen(ch);
	int last = -1, id = 0;
	nxt[0] = -1;
	for (int i = 1; i < n; i++) {
		int j;
		if (last < i) j = 0;
		else j = min(nxt[i - id], last - i) + 1;
		while (j < n && ch[j] == ch[i + j]) j++;
		nxt[i] = j - 1;
		if (i + nxt[i] > last) {
			id = i;
			last = i + nxt[i];
		}
	}
}

int Left[maxn], Right[maxn];

int main() {
	int t;
	scanf("%d", &t);
	while (t--) {
		scanf("%s", ch);
		n = strlen(ch);
		for (int i = 0; i < n; i++) {
			exkmp(ch + i);
			for (int j = i + 1; j < n; j++) 
				if (nxt[j - i] >= j - i - 1) {
					Left[2 * j - i]++;
				}
		}
		for (int i = 0; 2 * i < n; i++) swap(ch[i], ch[n - i - 1]);
		for (int i = 0; i < n; i++) {
			exkmp(ch + i);
			for (int j = i + 1; j < n; j++) 
				if (nxt[j - i] >= j - i - 1) {
					Right[n - (2 * j - i)]++;
				}
		}
		long long ans = 0;
		for (int i = 0; i < n; i++) ans += 1LL * Left[i] * Right[i];
		for (int i = 1; i <= n; i++) Left[i] = Right[i] = 0;
		printf("%lld\n", ans);
	}
} 
```
