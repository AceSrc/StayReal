
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
