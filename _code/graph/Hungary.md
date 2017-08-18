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
