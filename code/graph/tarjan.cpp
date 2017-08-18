
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
