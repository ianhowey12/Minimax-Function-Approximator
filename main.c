#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

typedef struct {
    double errorMin;
    double errorMax;
    double errorAv;
    double coef[100];
} approx;

/*
Find the best approximation in [A, B] with the given term form for the given function.
Stores the approximation data in D if it finds an approximation close to f
at frequency x-values spaced evenly in each interval between the points
where the approximation intersects f.
term is the function that takes an x-value and a vector to store the terms in.
N is the number of terms in the numerator.
M is the number of terms in the denominator.
points are the starting points at which the approximation will initially intersect f.
stepScale and stepMax between 0 and 1 determine the algorithm step size.
*/
bool approximate(approx* D, double f(double), double term(double, double*), int N, int M, double A, double B, int frequency, double* points, int iterations, double stepScale, double stepMax, bool print){
    if(N < 0 || N > 100) return 0;
    if(M < 0 || M > 100) return 0;
    int K = N + M;
    if(K > 100) return 0;
    if(A >= B) return 0;
    if(frequency < 2) return 0;
    if(stepScale < 0.0 || stepScale > 1.0) return 0;
    if(stepMax < 0.0 || stepMax > 1.0) return 0;
    if(iterations < 0) return 0;

    for(int i=0;i<K;i++){
        double x = points[i];
        if(x < A || x > B) return 0;
    }
    for(int i=1;i<K;i++){
        if(points[i-1] >= points[i]) return 0;
    }

    for(int T=0;T<iterations;T++){

        // Create the K x K+1 matrix.
        double v[100][101];
        for(int i=0;i<K;i++){
            double x = points[i];
            double y = f(x);

            // Compute the vector result of the terms.
            double r[100];
            term(x, r);

            // Set numerator elements as the computed term and denominator elements as y times that.
            for(int j=0;j<K;j++){
                double z = r[j];
                if(j >= N){
                    z *= y;
                }
                v[i][j] = z;
            }
            v[i][K] = y;
        }

        // Row reduce the matrix to solve for the coefficients.
        for(int i=0;i<K;i++){ // iterate over columns
            // Zero out the other elements in this column by adding a multiple of this row to other rows.
            for(int j=0;j<K;j++){ // iterate over rows
                if(i == j) continue;
                if(v[i][i] == 0.0){
                    printf("Iteration %i: interpolation problem could not be solved.\n", T);
                    return 0;
                }
                double ratio = v[j][i] / v[i][i];
                v[j][i] = 0.0;
                for(int a=i+1;a<=K;a++){ // iterate over row elements
                    v[j][a] -= v[i][a] * ratio;
                }
            }
        }

        // Normalize the pivot elements.
        for(int i=0;i<K;i++){ // iterate over rows
            double s = 1.0 / v[i][i];
            v[i][i] = 1.0;
            for(int j=i+1;j<=K;j++){ // iterate over row elements
                v[i][j] *= s;
            }
        }

        // Store the results (last column).
        for(int i=0;i<K;i++){
            if(i < N){
                D->coef[i] = v[i][K];
            }else{
                D->coef[i] = -v[i][K];
            }
        }
        
        if(print){
            printf("ITERATION %i\n\n", T);
            printf("PARAMETERS\n");
            for(int i=0;i<K;i++){
                printf("%.20lf  ", D->coef[i]);
            }
            printf("\n\nERRORS\n");
        }
        
        // Find the intervals between the test points.
        // Find the maximum error in the intervals between the test points.
        double e[101];
        double ex[101];
        for(int c=0;c<=K;c++){
            e[c] = 0.0;
            ex[c] = 0.0;

            double x0 = A, x1 = B;
            if(c > 0) x0 = points[c-1];
            if(c < K) x1 = points[c];

            for(int i=0;i<frequency;i++){
                
                double x = x0 + (x1 - x0) * i / (double)(frequency - 1);
                if(i == frequency - 1) x = x1;

                // Compute the vector result of the terms.
                double r[100];
                term(x, r);

                double numer = 0.0;
                double denom = 1.0;
                for(int j=0;j<K;j++){
                    if(j < N){
                        numer += D->coef[j] * r[j];
                    }else{
                        denom += D->coef[j] * r[j];
                    }
                }

                double error = (numer / denom) - f(x);
                if(error < 0.0) error = -error;
                if(error > e[c]){
                    e[c] = error;
                    ex[c] = x;
                }
            }

            if(print){
                printf("%.20lf at %.5lf (interval %.5lf)\n", e[c], ex[c], x1 - x0);
            }
        }

        // Find the average, min, and max of all the errors.
        D->errorMin = 1e300;
        D->errorMax = -1e300;
        D->errorAv = 0.0;
        for(int c=0;c<=K;c++){
            double x = e[c];
            if(x > D->errorMax){
                D->errorMax = x;
            }
            if(x < D->errorMin){
                D->errorMin = x;
            }
            D->errorAv += e[c];
        }
        D->errorAv /= (double)(K + 1);

        // Find the maximum deviation.
        double md = 0.0;
        for(int c=0;c<=K;c++){
            double x = e[c] - D->errorAv;
            if(x < 0.0) x = -x;
            if(x > md) md = x;
        }

        // Find the normalized deviations.
        double d[101];
        for(int c=0;c<=K;c++){
            d[c] = (e[c] - D->errorAv) / md;
        }

        // Calculate the step size.
        double s = stepScale * md / D->errorAv;
        if(s > stepMax) s = stepMax;

        // Find the scaled interval lengths.
        // Adding 1 to a normalized deviation multiplies that unnormalized interval length by 1.0 - step.
        // Therefore, intervals with large errors get shrunk (since 1.0 - step < 1).
        double l[101];
        for(int c=0;c<=K;c++){
            double x0 = A, x1 = B;
            if(c > 0) x0 = points[c-1];
            if(c < K) x1 = points[c];

            double curr = e[c];
            l[c] = pow(1 - s, d[c]) * (x1 - x0);
        }

        // Sum the lengths.
        double total = 0.0;
        for(int c=0;c<=K;c++){
            total += l[c];
        }

        // Make sure every length is sufficiently large.
        for(int c=0;c<=K;c++){
            l[c] /= total;
            double epsilon = 1e-15;
            if(l[c] < epsilon){
                l[c] = epsilon;
            }
        }

        // Take the prefix sum of lengths.
        double sum = 0.0;
        for(int c=0;c<=K;c++){
            sum += l[c];
            l[c] = sum;
        }

        // Update the test points.
        for(int c=0;c<K;c++){
            points[c] = A + (B - A) * l[c] / l[K];
        }

        if(print){
            printf("POINTS: ");
            for(int c=0;c<K;c++){
                printf("%.20lf ", points[c]);
            }
            printf("\n\n");
        }
    }

    if(print){
        printf("PARAMETERS\n");
        for(int i=0;i<K;i++){
            printf("%.20lf  ", D->coef[i]);
        }
        printf("\n\n");

        printf("%.20lf %.20lf\n", D->errorMin, D->errorMax);
    }

    return 1;
}
