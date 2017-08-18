- math
  - [fft](#fft)
- string
  - [sam](#sam)
  - [z-algorithm](#zalgorithm)
  - [SuffixArray](#suffixarray)
  - [ac_automation](#ac_automation)
- number_theory
  - [FWT](#fwt)
  - [ntt](#ntt)
- graph
  - [Hungary](#hungary)
  - [costFlow](#costflow)
  - [tarjan](#tarjan)
  - [maxFlow](#maxflow)
  - [cutPoint](#cutpoint)
  - [bridge](#bridge)
  - [dmst](#dmst)
- data_structure
  - [Splay](#splay)
  - [virtualTree](#virtualtree)


# fft

```cpp

const int N = 1000000;
const int maxn = 1000000 + 2;
const double pi = 3.1415926535897932384626433832795;

struct Point {
	double x, y;
	Point (double a = 0, double b = 0) : x(a), y(b) {
	}
} x1[N / 2], x2[N / 2];

char a[N/2],b[N/2];  
int sum[N]; 

Point operator + (Point &a, Point &b) {
	return Point(a.x + b.x, a.y + b.y);
}

Point operator - (Point &a, Point &b) {
	return Point(a.x - b.x, a.y - b.y);
}

Point operator * (Point &a, Point &b) {
	return Point(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

void trans(Point *p, int bit) {
	int n = (1 << bit), i, j = (1 << bit) >> 1;
	for (i = 1; i < n - 1; i++) {
		if (i < j) swap(p[i], p[j]);
		int t = (1 << bit) >> 1;
		while (j >= t) {
			j -= t;
			t >>= 1;
		}
		j += t;
	}
}

void fft(Point *y, int bit, int on) {
	trans(y, bit);
	
	int n = (1 << bit);
	for (int s = 0; s < bit; s++) {
		int m = (1 << s) << 1, p = (1 << s);
		Point w = Point(cos(on * 2 * pi / m), sin(on * 2 * pi / m));
		for (int i = 0; i < n; i += m) {
			Point cur = Point(1, 0);
			for (int j = 0; j < p; j++) {
				Point t = y[i + j];
				Point s = y[i + j + p] * cur;
				y[i + j] = t + s;
				y[i + j + p] = t - s;
				cur = cur * w;
			}
		}
		//for (int i = 0; i < n; i++) printf("%.6lf\n", y[i].x);
		//printf("\n");
	}
	if (on == -1) 
		for (int i = 0; i < n; i++) y[i].x /= n;
	//printf("--------------\n");
}

int ans[maxn];

int main() {

	int n, m;
	scanf("%d%d", &n, &m);
	n++;
	m++;
	for (int i = 0; i < n; i++)
		scanf("%lf", &x1[i].x);
	for (int i = 0; i < m; i++)
		scanf("%lf", &x2[i].x);
	int bit = 0, least = max(n, m);
	while ((1 << bit) < (least << 1)) bit++;
	fft(x1, bit, 1);
	fft(x2, bit, 1);
	for (int i = (1 << bit) - 1; i > -1; i--) x1[i] = x1[i] * x2[i];
	fft(x1, bit, -1);
	n = n + m - 1;
	for (int i = 0; i < n; i++) ans[i] = x1[i].x + 0.5;
	for (int i = 0; i < n - 1; i++) printf("%d ", ans[i]);
	printf("%d\n", ans[n - 1]);
}  

```


# sam

```cpp
#include <stdio.h>
#include <string.h>

const int maxn = 1000000;
const char c = 26;

struct Tnode {
	Tnode *nxt[c], *par;
	int val;
} nodes[maxn], *root, *last;

int N;
Tnode *newNode(int x) {
	nodes[N].val = x;
	N++;
	return nodes + N - 1;
}

void insert(char s) {
	Tnode *p = last;
	Tnode *np = newNode(p->val + 1);
	while (p && p->nxt[s] == 0) p->nxt[s] = np, p = p->par;
	if (!p) np->par = root;
	else {
		Tnode *q = p->nxt[s];
		if (p->val + 1 == q->val) {
			np->par = q;
		} 
		else {
			Tnode *nq = newNode(p->val + 1);
			memcpy(nq->nxt, q->nxt, sizeof(q->nxt));
			nq->par = q->par;
			q->par = nq;
			np->par = nq;
			while (p && p->nxt[s] == q)
				p->nxt[s] = nq, p = p->par; 
		}
	}
	last = np;
}

bool accept(char *s) {
	int n = strlen(s);
	Tnode *u = root;
	for (int i = 0; i < n; i++) {
		if (!u) return false;
		u = u->nxt[ s[i] - 'a' ];
	}
	if (!u) return false;
	return true;
}

char s[maxn];
int main() {
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	scanf("%s", s);
	int n = strlen(s);
	last = root = newNode(0);
	for (int i = 0; i < n; i++) insert(s[i] - 'a');
	scanf("%s", s);
	n = strlen(s);
	
	Tnode *u = root;
	int t = 0, rt = 0;
	for (int i = 0; i < n; i++) {
		char c = s[i] - 'a';
		if (u->nxt[c]) u = u->nxt[c], t++;
		else {
			while (u && !u->nxt[c]) u = u->par;
			if (!u) u = root, t = 0;
			else t = u->val + 1, u = u->nxt[c];
		}
		if (t > rt) rt = t;
	}
	printf("%d\n", rt);
	return 0;
}


```


# zalgorithm

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


# suffixarray

```cpp

void work() {
    int *x = ra, *y = rb, *t, m = 27, i, j, p;
    for (i = 0; i < n; i++) sum[ x[i] = r[i] ]++;
    for (i = 1; i < m; i++) sum[i] += sum[i - 1];
    for (i = n - 1; i > -1; i--) sa[ --sum[ x[i] ] ] = i;

    for (p = 0, j = 1; p < n; j <<= 1) {
        for (p = 0, i = n - j; i < n; i++) y[p++] = i;
        for (i = 0; i < n; i++)
            if (sa[i] >= j) y[p++] = sa[i] - j;

        for (i = 0; i < m; i++) sum[i] = 0;
        for (i = 0; i < n; i++) sum[ x[i] ]++;
        for (i = 1; i < m; i++) sum[i] += sum[i - 1];
        for (i = n - 1; i > -1; i--) sa[ --sum[ x[ y[i] ] ] ] = y[i];

        t = x; x = y; y = t;
        p = 1; x[ sa[0] ] = 0;
        for (i = 1; i < n; i++)
            if (cmp(y, sa[i - 1], sa[i], j)) x[ sa[i] ] = p - 1;
            else x[ sa[i] ] = p++;
        m = p;
    }
    n--;
    for (i = 0; i < n; i++) sa[i] = sa[i + 1];
}

void getHeight() {
    for (int i = 0; i < n; i++) rank[ sa[i] ] = i;
    int p = 0;
    for (int i = 0; i < n; i++) {
        if (p != 0) p--;
        if (rank[i] == 0) continue;
        int j = sa[ rank[i] - 1 ];
        while (r[i + p] == r[j + p]) p++;
        height[ rank[i] ] = p;
    }
}

```


# ac_automation

```cpp
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

```


# fwt

```cpp
xor:
	tf = (tf(A0)+tf(A1), tf(A0)-tf(A1))
	utf = (utf((A0+A1) / 2), utf((A0-A1)/ 2))
and:
	tf = (tf(A0)+tf(A1), tf(A1));
	utf = (utf(A0)-utf(A1), utf(A1));
or:
	tf = (tf(A0), tf(A1)+tf(A0))
	utf = (utf(A0), utf(A1)-utf(A0))

void FWT( ll X[], int l, int r, int v ) {
    if ( l == r ) return;
    int m = ( l + r ) >> 1;
    FWT( X, l, m, v ); FWT( X, m + 1, r, v );
    FOR ( i, 0, m - l ) {
        X[ l + i ] += X[ m + 1 + i ] * v;
    }
```


# ntt

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


# hungary

```cpp
#include <stdio.h>

const int maxt = 250001;
const int maxn = 501;

int n, m, t;
int boy[maxn], girl[maxn];
bool map[maxn][maxn], vis[maxn];

bool trace(int u, int v) {
    if (vis[v]) return false;
    vis[v] = true;
    if (!girl[v]) {
        boy[u] = v;
        girl[v] = u;
        return true;
    }
    for (int i = 1; i <= m; i++)
        if (map[ girl[v] ][i] && trace(girl[v], i)) {
            boy[u] = v;
            girl[v] = u;
            return true;
        }
    return false;
}

int main() {
    int a, b;
    scanf("%d%d%d", &n, &m, &t);
    while (t--) {
        scanf("%d%d", &a, &b);
        map[a][b] = true;
    }
    
    int ans = 0;
    for (int i = 1; i <= n; i++) {
        if (boy[i]) continue;
        for (int j = 1; j <= m; j++) vis[j] = false;
        for (int j = 1; j <= m; j++)
            if (map[i][j] && trace(i, j)) {
                ans++;
                break;
            }
    }
    printf("%d\n", ans);
    for (int i = 1; i <= n; i++) printf("%d ", boy[i]);
    printf("\n");
    return 0;
}

```


# costflow

```cpp

```


# tarjan

```cpp

void tarjan(int u, int pre) {
	dfn[u] = low[u] = ++times;
	state[u] = 1;
	for (int i = head[u]; i != -1; i = e[i].next) {
		int v = e[i].v;
		if (v == pre) continue;
		if (!state[v]) {
			tarjan(v, u);
			low[u] = min(low[u], low[v]);
		}
		else if (state[v] == 1) low[u] = min(low[u], dfn[v]);
	}
	if (dfn[pre] <= low[u]) {
		if (pre == 1) {
			rootCnt++;
		}
		else cut[pre] = true;
	}
}

```


# maxflow

```cpp
namespace Network {
    const int maxn = 600002;
    const int maxm = 600002;
    const int inf = 2147483647;

    struct Tedge {
        int v, next, f;
    } e[maxm * 2];

    int N;
    int head[maxn];

    void insert(int x, int y, int c) {
        e[N].v = y; e[N].next = head[x]; e[N].f = c;
        head[x] = N; N++;
        if (N & 1) insert(y, x, 0);
    }

    int st, ed;
    int o[maxn], dis[maxn];

    void init(int S, int T) {
        st = S; ed = T;
        for (int i = S; i <= T; i++) head[i] = -1;
        N = 0;
    }

    bool bfs() {
        for (int i = 1; i <= ed; i++) dis[i] = -1;
        o[0] = st; dis[st] = 0;
        for (int h = 0, t = 0; h <= t; h++) {
            int u = o[h];
            if (u == ed) return true;
            for (int i = head[u]; i != -1; i = e[i].next)
                if (e[i].f > 0 && dis[ e[i].v ] == -1) {
                    o[++t] = e[i].v;
                    dis[ e[i].v ] = dis[u] + 1;
                }
        }
        return false;
    }

    int dfs(int u, int limit) {
        if (u == ed) return limit;
        int rt = 0;
        for (int i = head[u]; i != -1; i = e[i].next) 
            if (dis[ e[i].v ] == dis[u] + 1 && e[i].f > 0) {
                int tmp = dfs(e[i].v, min(e[i].f, limit - rt));
                e[i].f -= tmp;
                e[i ^ 1].f += tmp;
                rt += tmp;
                if (rt == limit) return rt;
            }	
        dis[u] = -1;
        return rt;
    }

    int maxFlow() {
        int ans = 0;
        while (bfs()) ans += dfs(st, inf);
        return ans;
    }
};


```


# cutpoint

```cpp
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
```


# bridge

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


# dmst

```cpp
#include <vector>
#include <cassert>
#include <cstdio>
#include <algorithm>
using namespace std;

#define FOR(i, a, b) for (int i = (a); i < (b); ++i)
#define REP(i, n) FOR(i, 0, n)

typedef long long llint;

struct Edge {
  int x, y, w;
};

bool operator < (const Edge &a, const Edge &b) {
  return a.x < b.x || (a.x == b.x && a.y < b.y);
}

int dmst(int N, vector <Edge> &E, int root) {
  const int inf = -1e9;
  vector<int> cost(N), back(N), label(N), bio(N);
  int ret = 0;
  while (1) {
    REP(i, N) cost[i] = inf;
    for (auto e:  E) {
      if (e.x == e.y) continue;
      if (e.w > cost[e.y]) cost[e.y] = e.w, back[e.y] = e.x;
    }
    cost[root] = 0;
    REP(i, N) if (cost[i] == inf) return -1;
    REP(i, N) ret += cost[i];
    int K = 0;
    REP(i, N) label[i] = -1;
    REP(i, N) bio[i] = -1;
    REP(i, N) {
      int x = i;
      for (; x != root && bio[x] == -1; x = back[x])
        bio[x] = i;
      if (x != root && bio[x] == i) {
        for (; label[x] == -1; x = back[x])
          label[x] = K;
        ++K;
      }
    }
    if (K == 0) break;
    REP(i, N) if (label[i] == -1)
      label[i] = K++;

    for (auto & e:  E) {
      int xx = label[e.x];
      int yy = label[e.y];
      if (xx != yy)
        e.w -= cost[e.y];
      e.x = xx;
      e.y = yy;
    }
    root = label[root];
    N = K;
  }
  return ret;
}

void check_connectivity(int N, vector<Edge> &E) {
  vector<bool> vis(N);
  vis[0] = true;
  REP(i, N) {
    for (auto e: E) {
      if (vis[e.x]) vis[e.y] = true;
    }
  }
  REP(i, N) {
    assert(vis[i]);
  }
}

void check_dup(vector<Edge> &E) {
  sort(E.begin(), E.end());
  REP(i, E.size() - 1) {
    Edge &a = E[i];
    Edge &b = E[i + 1];
    assert(a.x != b.x || a.y != b.y);
  }
}

int main() {
  int tn, n, m;
  scanf("%d", &tn);
  REP(ti, tn) {
    scanf("%d%d", &n, &m);
    assert(1 <= n && n <= 1000 && 1 <= m && m <= 10000);
    vector<Edge> edges;
    REP(i, m) {
      Edge e;
      scanf("%d%d%d", &e.x, &e.y, &e.w);
      assert(1 <= e.x && e.x <= n);
      assert(1 <= e.y && e.y <= n);
      assert(1 <= e.w && e.w <= 100);
      assert(e.y != 1);
      e.x -= 1; e.y -= 1;
      e.w *= n;
      if (e.y == n - 1) {
        e.w += (n - 1 - e.x);
      }
      edges.push_back(e);
    }
    check_connectivity(n, edges);
    check_dup(edges);
    int ans = dmst(n, edges, 0);
    int cost = ans / n;
    int idx = n - 1 - (ans % n);
    printf("%d %d\n", cost, idx + 1);
  }
  return 0;
}

```


# splay

```cpp
#include <stdio.h>

const int maxn = 1000002;

struct Tnode {
	Tnode *ch[2], *par;
	int dir, best, size, num;
} nodes[maxn], *root;

int N;
Tnode *newNode(int no) {
	// U have to update par.
	nodes[N].best = nodes[N].num = no;
	nodes[N].ch[0] = nodes[N].ch[1] = 0;
	nodes[N].size = nodes[N].dir = 1;
	N++;
	return nodes + N - 1;
}

void up(Tnode *u) {
	u->size = 1; u->best = u->num;
	if (u->ch[0]) {
		u->size += u->ch[0]->size;
		if (u->ch[0]->best > u->best) u->best = u->ch[0]->best;
	}
	if (u->ch[1]) {
		u->size += u->ch[1]->size;
		if (u->ch[1]->best > u->best) u->best = u->ch[1]->best;
	}
}

void rotate(Tnode *u) {
	Tnode *par = u->par;
	char dir = u->dir;
	u->par = par->par; u->dir = par->dir;
	if (u->par) u->par->ch[ par->dir ] = u;
	par->ch[dir] = u->ch[dir ^ 1];
	if (u->ch[dir ^ 1]) par->ch[dir]->par = par, par->ch[dir]->dir = dir;
	u->ch[dir ^ 1] = par;
	par->par = u;
	par->dir = dir ^ 1;
	up(par);
}

void splay(Tnode *u, Tnode *tar) {
	while (u->par != tar) {
		if (u->par->par != tar) 
			if (u->dir == u->par->dir) rotate(u->par);
			else rotate(u);
		rotate(u);
	}
	up(u); up(tar); // tar is not updated.
}

Tnode *getNode(int k) {
	Tnode *u = root->ch[1];
	int c;
	while (k) {
		if (!u->ch[0]) c = 0;
		else c = u->ch[0]->size;
		if (c + 1 == k) return u;
		if (c < k) {
			k -= c + 1;
			u = u->ch[1];
		}
		else u = u->ch[0];
	}
}

void merge(Tnode *a, Tnode *b) {
	splay(a, root);
	splay(b, a);
}

void merge(int x, int y) {
	merge(getNode(x), getNode(y + 2));
}

void build(int n, int *d) {
	N = 0;
	root = newNode(0); 
	Tnode *LeftGuard = newNode(0), *RightGuard = newNode(0), *cur, *tmp;
	LeftGuard->par = root; cur = LeftGuard;
	for (int i = 1; i <= n; i++) {
		tmp = newNode(d[i]);
		tmp->par = cur; cur = tmp;
	}
	RightGuard->par = tmp;
	splay(RightGuard, root);
}

int query(int a, int b) {
	merge(a, b);
	return root->ch[1]->ch[1]->ch[0]->best;
}

int alter(int a, int b) {
	splay(getNode(a + 1), root);
	root->ch[1]->num = b;
	splay(root->ch[1], root);
}

int d[maxn];
int main() {
	freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
	int n, m, a, b;
	while (scanf("%d%d", &n, &m) != EOF) {
		for (int i = 1; i <= n; i++) scanf("%d", d + i);
		build(n, d);
		char ch;
		while (m--) {
			for (ch = getchar(); ch != 'Q' && ch != 'U'; ch = getchar());
			scanf("%d%d", &a, &b);
			if (ch == 'Q') printf("%d\n", query(a, b));
			else alter(a, b);
		}
	}	
	return 0;
}

```


# virtualtree

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
