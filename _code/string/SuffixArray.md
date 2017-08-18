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
