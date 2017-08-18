#include <stdio.h>
#include <string.h>
#include <algorithm>

const int maxn = 70002;
const int c = 128;

int root;
int nxt[maxn][c], fail[maxn], s[maxn];
int N;
int newNode() {
	N++;
	return N;
} 

void insert(char *_s, int id) {
	int n = strlen(_s);
	int u = root;
	for (int i = 0; i < n; i++) {
		char c = _s[i];
		if (!nxt[u][c]) nxt[u][c] = newNode();
		u = nxt[u][c];
	}
	s[u] = id;
}

int o[maxn];
void build() {
	int h = 0, t = 0;
	for (int i = 0; i < c; i++)
		if (nxt[0][i]) {
			o[t++] = nxt[0][i];
			fail[ nxt[0][i] ] = root;
		}
	while (h < t) {
		int u = o[h++];
		for (int i = 0; i < c; i++) 
			if (nxt[u][i]) {
				fail[ nxt[u][i] ] = nxt[ fail[u] ][i];
				o[t++] = nxt[u][i];
			}
			else nxt[u][i] = nxt[ fail[u] ][i];
	}
}

char _s[maxn];
int list[maxn];

bool work(char *_s, int id) {
	int n = strlen(_s), t = 0, cnt = 0;
	int u = root, p;
	int rt = 0;
	for (int i = 0; i < n; i++) {
		char c = _s[i];
		u = nxt[u][c]; p = u;
		while (p) {
			if (s[p]) list[cnt++] = s[p];
			p = fail[p];
		}
	}
	if (cnt == 0) return false;
	printf("web %d: ", id);
	std :: sort(list, list + cnt);
	printf("%d", list[0]);
	for (int i = 1; i < cnt; i++) 
		if (list[i] != list[i - 1]) printf(" %d", list[i]);
	putchar('\n');
	return true;
}

int main() {
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	int n, cnt = 0;
	scanf("%d", &n);
	for (int i = 1; i <= n; i++) {
		scanf("%s", _s);
		insert(_s, i);
	}
	build();
	scanf("%d", &n);
	for (int i = 1; i <= n; i++) {
		scanf("%s", _s);
		if (work(_s, i)) cnt++;;
	}
	printf("total: %d\n", cnt);
	return 0;
} 
