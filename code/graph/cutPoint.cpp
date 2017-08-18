#include <fstream>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <string>
#include <stack>
using namespace std;

const int maxn = 1000 + 2;

struct Tedge {
	int next, v;
} e[maxn * maxn];

int head[maxn], N;

void Insert(int x, int y) {
	e[N].v = y;
	e[N].next = head[x];
	head[x] = N;
	N++;
}

char ch;

int get() {
	ch = ' ';
	int rt = 0;
	while (ch < '0' || ch > '9') ch = getchar();
	for ( ; ch >= '0' && ch <= '9'; ch = getchar()) rt = rt * 10 + ch - '0';
	return rt;
}

int n, t;
int dfn[maxn], low[maxn], state[maxn], anc[maxn];
bool isCut[maxn];

bool Init() {
	n = get();
	if (n == 0) return false;
	for (int i = 0; i < n; i++) head[i] = -1;
	while (true) {
		int u = get();
		if (u == 0) break;
		u--;
		while (true) {
			int v = get();
			v--;
			Insert(u, v);
			Insert(v, u);
			if (ch == '\n') break;
		}
	}
	return true;
}

void Tarjan(int u, int par) {
	dfn[u] = t;
	low[u] = t;
	state[u] = 1;
	anc[u] = par;
	t++;
	int cnt = 0;
	for (int i = head[u]; i != -1; i = e[i].next) {
		int v = e[i].v;
		if (state[v] == 0) {
			Tarjan(v, u);
			low[u] = min(low[u], low[v]);
			cnt++;
		}
		else if (state[v] == 1) low[u] = min(low[u], dfn[v]);
	}
	if (u == 0 && cnt > 1) {
		isCut[u] = true;
		return ;
	}
	if (u != 0) {
		for (int i = head[u]; i != -1; i = e[i].next) {
			int v = e[i].v;
			if (anc[v] != u) continue;
			if (low[v] >= dfn[u]) {
				isCut[u] = true;
				return ;
			}
		}
	}
}

void Solve() {
	for (int i = 0; i < n; i++) state[i] = 0;
	for (int i = 0; i < n; i++) isCut[i] = false;
	Tarjan(0, -1);
	int s = 0;
	for (int i = 0; i < n; i++)
		if (isCut[i]) s++;
	cout << s << endl;
}

int main() {
	while (Init()) 
		Solve();
	return 0;
}