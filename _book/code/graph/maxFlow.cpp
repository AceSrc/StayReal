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

