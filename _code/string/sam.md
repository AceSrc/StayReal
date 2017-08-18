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
