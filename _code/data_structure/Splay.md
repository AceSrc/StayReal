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
