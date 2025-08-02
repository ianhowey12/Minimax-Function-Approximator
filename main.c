#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* inverseZTable;
int inverseZTableSize;

double* extremeTable;
int extremeTableSize;
double extremeTableThreshold;


double zToProb(double z){
    return 0.5 * (1.0 + erf(z / sqrt(2.0)));
}

// slow but perfectly accurate
double findZGivenProb(double p, double l, double r){
    double m = (l + r) / 2;
    if(r - l <= 0.00000000001) return m;
    if(zToProb(m) > p) return findZGivenProb(p, l, m);
    return findZGivenProb(p, m, r);
}

void generateInverseZTable(int numValuesInTable, int numValuesInExtremeTable, double threshold){
    if(numValuesInTable < 10 || numValuesInTable > 100000000){
        printf("The number of values to generate the inverse z-table must be between 10 and 100000000.\n\n");
        exit(1);
    }

    inverseZTableSize = numValuesInTable;
    inverseZTable = (double*)malloc(sizeof(double) * (numValuesInTable + 1));

    // use 2 values before last to find last value
    double upperBound = 2 * findZGivenProb((double)(numValuesInTable - 1) / (double)(numValuesInTable), -7.5, 7.5) - findZGivenProb((double)(numValuesInTable - 2) / (double)(numValuesInTable), -7.5, 7.5);

    for(int i = 0; i <= numValuesInTable; i++){
        inverseZTable[i] = findZGivenProb((double)(i) / (double)(numValuesInTable), -upperBound, upperBound);
    }

    extremeTableSize = numValuesInExtremeTable;
    extremeTableThreshold = threshold;
    extremeTable = (double*)malloc(sizeof(double) * (numValuesInExtremeTable + 1));

    // use 2 values before last to find last value
    upperBound = 2 * findZGivenProb(threshold + (1.0 - threshold) * (double)(numValuesInExtremeTable - 1) / (double)(numValuesInExtremeTable), -7.5, 7.5) - findZGivenProb(threshold + (1.0 - threshold) * (double)(numValuesInExtremeTable - 2) / (double)(numValuesInExtremeTable), -7.5, 7.5);

    for(int i = 0; i <= numValuesInExtremeTable; i++){
        extremeTable[i] = findZGivenProb(threshold + (1.0 - threshold) * (double)(i) / (double)(numValuesInExtremeTable), -upperBound, upperBound);
    }
}

// fast and less accurate
double findZFromTable(double p){
    if(p >= extremeTableThreshold){
        double i = (p - extremeTableThreshold) / (1.0 - extremeTableThreshold) * (double)extremeTableSize;
        double l = extremeTable[(int)i];
        double r = extremeTable[extremeTableSize];
        if(i < extremeTableSize) r = extremeTable[(int)i + 1];
        // linear approximation between l and r
        return l + (r - l) * (i - (int)i);
    }
    if(p <= 1.0 - extremeTableThreshold){
        double i = (1.0 - extremeTableThreshold - p) / (1.0 - extremeTableThreshold) * (double)extremeTableSize;
        double l = extremeTable[(int)i];
        double r = extremeTable[extremeTableSize];
        if(i < extremeTableSize) r = extremeTable[(int)i + 1];
        // linear approximation between l and r
        return -(l + (r - l) * (i - (int)i));
    }
    double i = p * (double)inverseZTableSize;
    double l = inverseZTable[(int)i];
    double r = inverseZTable[inverseZTableSize];
    if(i < inverseZTableSize) r = inverseZTable[(int)i + 1];
    // linear approximation between l and r
    return l + (r - l) * (i - (int)i);
}

void printInverseZTable(){
    printf("Inverse Z-Table (size %i):\n", inverseZTableSize);
    for(int i=0;i<=inverseZTableSize;i++){
        printf("%-20.14f%-20.14lf\n", (double)i / (double)inverseZTableSize, inverseZTable[i]);
    }
}

void printExtremeTable(){
    printf("\nInverse Z-Table Extreme Values (size %i, threshold %.14f):\n", extremeTableSize, extremeTableThreshold);
    for(int i=0;i<=extremeTableSize;i++){
        printf("%-20.14f%-20.14lf\n", ((double)i / (double)extremeTableSize) * (1.0 - extremeTableThreshold) + extremeTableThreshold, extremeTable[i]);
    }
    printf("\n");
}

void testAccuracy(int numTrials, char printAllResults){
    if(numTrials < 1 || numTrials > 1000000){
        printf("The number of trials in the accuracy test must be between 1 and 1000000.\n\n");
        exit(1);
    }

    double avDev = 0.0;
    double var = 0.0;
    double sd = 0.0;
    if(printAllResults) printf("%-20s%-20s%-20s%-20s\n", "Probability", "Real Z-Score", "Approx. Z-Score", "Deviation");
    for(int i=0;i<numTrials;i++){
        double p = ((double)rand() * 32768.0 + (double)rand()) / (32768.0 * 32768.0);
        double a = findZGivenProb(p, -7.5, 7.5);
        double b = findZFromTable(p);

        double dev = fabs(a - b);
        avDev += dev;
        var += dev * dev;
        
        if(printAllResults) printf("%-20.14lf%-20.14lf%-20.14lf%-20.14lf\n", p, a, b, dev);
    }
    avDev /= (double)numTrials;
    var /= (double)numTrials;
    sd = sqrt(var);
    printf("Average Deviation: %.14f, Average Squared Deviation (Variance): %.14f, Standard Deviation: %.14f, Number of Trials: %i\n\n", avDev, var, sd, numTrials);
}

void testExtremeValues(){
    double extremes[] = {0.0, 0.00000001, 0.0001, 1.0 - extremeTableThreshold - 0.00000001, 1.0 - extremeTableThreshold, 0.5, extremeTableThreshold - 0.00000001, extremeTableThreshold, 0.9999, 0.99999999, 1.0};

    for(int i=0;i<11;i++){
        double p = extremes[i];
        double a = findZGivenProb(p, -7.5, 7.5);
        double b = findZFromTable(p);

        double dev = fabs(a - b);
        printf("%-20.14lf%-20.14lf%-20.14lf%-20.14lf\n", p, a, b, dev);
    }
}

int main(void){

    generateInverseZTable(100000, 100000, 0.9);

    //printInverseZTable();
    //printExtremeTable();

    testAccuracy(100000, 0);
    testExtremeValues();

    return 0;
}
