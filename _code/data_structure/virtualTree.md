```cpp
#include <fstream>
#include <iostream>
#include <algorithm>
#include <string.h>
#include <math.h>
using namespace std;

const int maxn = 1000000 + 2;
const int inf = 2000000001;

struct Tedge {
	int v, next;
} e[maxn << 1];

int N;
int head[maxn];

void insert(int x, int y) {
	e[N].v = y;
	e[N].next = head[x];
	head[x] = N;
	N++;
}

int n, times;
int p[maxn], o[maxn], dfn[maxn], anc[maxn], id[maxn], dep[maxn];

char ch;
void get(int &n) {
	n = 0;
	for (ch = getchar(); ch < '0' || ch > '9'; ch = getchar()) ;
	for ( ; ch >= '0' && ch <= '9'; ch = getchar()) n = n * 10 + ch - '0';
}

void dfs() {
	for (int i = 1; i <= n; i++) p[i] = head[i];
	int t = 0;
	o[t] = 1;
	dfn[1] = ++times;
	id[times] = 1;
	while (t > -1) {
		int u = o[t];
		if (p[u] == -1) {
			t--;
			continue;
		}
		int v = e[ p[u] ].v;
		p[u] = e[ p[u] ].next;
		if (v == anc[u]) continue;
		dfn[v] = ++times;
		id[times] = v;
		anc[v] = u;
		dep[v] = dep[u] + 1;
		o[++t] = v;		
	}
}

int cnt[maxn], wei[maxn], pathroot[maxn];

void buildLca() {
	for (int i = n; i; i--) {
		int v = id[i];
		cnt[v]++;
		cnt[ anc[v] ] += cnt[v];
		if (!wei[ anc[v] ] || cnt[ wei[ anc[v] ] ] < cnt[v]) wei[ anc[v] ] = v;
	}
	pathroot[1] = 1;
	for (int i = 2; i <= n; i++)
		if (wei[ anc[ id[i] ] ] == id[i]) pathroot[ id[i] ] = pathroot[ anc[ id[i] ] ];
		else pathroot[ id[i] ] = id[i];
}

int Lca(int x, int y) {
	while (pathroot[x] != pathroot[y]) {
		if (dep[ pathroot[x] ] > dep[ pathroot[y] ]) x = anc[ pathroot[x] ];
		else y = anc[ pathroot[y] ];
	}
	if (dep[x] > dep[y]) return y;
	return x;
}

void Init() {
	int a, b;
	get(n);
	for (int i = 1; i <= n; i++) head[i] = -1;
	for (int i = 1; i < n; i++) {
		get(a);
		get(b);
		insert(a, b);
		insert(b, a);
	}
	dfs();
	buildLca();
}

bool cmp(int a, int b) {
	return dfn[a] < dfn[b];
}

int size, root;
int list[maxn], vir[maxn], in[maxn];
bool cut[maxn];

void virTree(int n) {
	sort(list, list + n, cmp);
	int t = 0;
	size = 0;
	vir[ list[0] ] = 0;
	o[t] = list[0];
	in[size++] = list[0];
	for (int i = 1; i < n; i++) {
		int u = list[i];
		int lca = Lca(o[t], u), ori = t;
		while (t > -1 && dep[ o[t] ] > dep[lca]) t--;
		if (t > -1 && o[t] == lca) {
			o[++t] = u;
			vir[u] = lca;
		}
		else if (t > -1 && o[t] != lca) {
			vir[lca] = o[t];
			vir[u] = lca;
			if (ori != t) vir[ o[t + 1] ] = lca;
			o[++t] = lca;
			o[++t] = u;
			in[size++] = lca;
		}
		else if (t == -1) {
			vir[ o[0] ] = lca;
			vir[u] = lca;
			vir[lca] = 0;
			o[++t] = lca;
			o[++t] = u;
			in[size++] = lca;
		}
		in[size++] = u;
	}
	N = 0;
	for (int i = 0; i < size; i++) head[ in[i] ] = -1;
	for (int i = 0; i < size; i++) 
		if (vir[ in[i] ]) insert(vir[ in[i] ], in[i]);
	root = o[0];
}

long long depCnt[maxn];
int best[maxn], least[maxn];

void dp() {
	o[0] = root;
	for (int h = 0, t = 0; h <= t; h++) {
		int u = o[h];
		for (int i = head[u]; i != -1; i = e[i].next) o[++t] = e[i].v;
	}
	long long ans = 0;
	int far = -inf, near = inf;
	for (int i = size - 1; i > -1; i--) {
		int u = o[i];
		if (!cut[u]) depCnt[u] = 0, cnt[u] = 0, best[u] = -inf, least[u] = inf;
		else depCnt[u] = dep[u], cnt[u] = 1, best[u] = 0, least[u] = 0;
		for (int j = head[u]; j != -1; j = e[j].next) {
			int v = e[j].v;
			ans += 1LL * cnt[v] * (depCnt[u] - 1LL * cnt[u] * dep[u]) + 1LL * cnt[u] * (depCnt[v] - 1LL * cnt[v] * dep[u]);
			cnt[u] += cnt[v];
			depCnt[u] += depCnt[v];
			
			if (best[u] != -inf && best[v] != -inf) far = max(far, best[u] + best[v] + dep[v] - dep[u]);
			if (least[u] != inf && least[v] != inf) near = min(near, least[u] + least[v] + dep[v] - dep[u]);
			if (best[v] != -inf) best[u] = max(best[u], best[v] + dep[v] - dep[u]);
			if (least[v] != inf) least[u] = min(least[u], least[v] + dep[v] - dep[u]);
		}
	}
	if (near == inf) near = 0;
	if (far == -inf) far = 0;
	printf("%lld %d %d\n", ans, near, far);
}

void Solve() {
	int m, k;
	get(m);
	while (m--) {
		get(k);
		for (int i = 0; i < k; i++) {
			get(list[i]);
			cut[ list[i] ] = true;
		}
		virTree(k);
		dp();
		for (int i = 0; i < k; i++) cut[ list[i] ] = false;
	}
}

int main() {
	Init();
	Solve();

	return 0;
}

```
