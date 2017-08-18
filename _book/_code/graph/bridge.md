```cpp
#include <fstream>
#include <iostream>
#include <algorithm>
using namespace std;

const int maxm = 1000 + 2;
const int maxn = 1000 + 2;

struct Tedge {
	int v, next;
} e[maxm * 2];

int N;
int head[maxn];

void insert(int x, int y) {
	e[N].v = y;
	e[N].next = head[x];
	head[x] = N;
	N++;
}

int n, r;

void Init() {
	int x, y;
	scanf("%d%d", &n, &r);
	for (int i = 0; i < n; i++) head[i] = -1;
	for (int i = 0; i < r; i++) {
		scanf("%d%d", &x, &y);
		x--;
		y--;
		insert(x, y);
		insert(y, x);
	}
}

int top, times, flag;
int stack[maxn], dfn[maxn], low[maxn], inx[maxn];
bool vis[maxn];

void dfs(int u, int pre) {


	stack[++top] = u;
	dfn[u] = low[u] = times++;
	vis[u] = true;
	for (int i = head[u]; i != -1; i = e[i].next) {
		int v = e[i].v;
		if (pre == v) continue;
		if (vis[v]) low[u] = min(low[u], dfn[v]);
		else {
			dfs(v, u);
			low[u] = min(low[u], low[v]);
			
			if (low[v] > dfn[u]) {
				flag++;
				while (true) {
					inx[ stack[top] ] = flag;
					top--;
					if ( stack[top + 1] == v ) break;
				}
			}
		}
	}
}

bool map[maxn][maxn];

void Solve() {
	for (int i = 0; i < n; i++)
		if (!vis[i]) {
			dfs(i, -1);
			if (top > -1) flag++;
			while (top > -1) inx[ stack[top--] ] = flag;
		}
	for (int i = 0; i < n; i++)
		for (int j = head[i]; j != -1; j = e[j].next) {
			int v = e[j].v;
			if (inx[v] != inx[i]) map[ inx[i] ][ inx[v] ] = true;
		}
	if (flag == 1) {
		printf("0\n");
		return ;
	}
	
	int c = 0;
	for (int i = 1; i <= flag; i++) {
		int cur = 0;
		for (int j = 1; j <= flag; j++)
			cur += map[i][j];
		if (cur == 1) c++;
	}
	
	cout << (c + 1) / 2 << endl;
}

int main() {

	Init();
	Solve();
	return 0;
}
```
