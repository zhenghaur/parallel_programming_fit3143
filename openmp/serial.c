#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

/* max false positive rate set to 5% */
#define MAX_FP 0.05

/* max word length to be read from file set to 100 */
#define MAX_WORD_LENGTH 100

uint32_t APHash(char *word, int wordLength);

uint32_t JenkinsHash(char *word, int wordLength);

int readFile(const char *fileName, int fileLength, char ***pppUniqueWordList);

void readQueryFile(const char *fileName, int fileLength, char ***pppQueryList, int **ppQueryTruthList);

bool isWordInList(char *word, char **wordList, int wordListLength);

void insert(bool *bloomfilter, int filterSize, char *word, int wordLength);

int query(bool *bloomfilter, int filterSize, char *word, int wordlength);

int main()
{
    struct timespec start, end, startComp, endComp;

    // get start time for entire code
    clock_gettime(CLOCK_MONOTONIC, &startComp);

    // initializing all variables
    double time_taken;

    char *fileNames[3] = {"datasets/MOBY_DICK.txt", "datasets/LITTLE_WOMEN.txt", "datasets/SHAKESPEARE.txt"};
    int fileLengths[3] = {215724, 195467, 965465};
    char *queryName = "datasets/query.txt";
    int queryLength = 91640;

    char **ppUniqueWordList[3] = {0};
    int uniqueWordListLength[3] = {0};

    char **ppQueryWordList = NULL;
    int *pQueryWordTruthList = NULL;

    int *results = (int *)malloc(queryLength * sizeof(int));

    int truePositive = 0;
    int falsePositive = 0;
    int trueNegative = 0;
    int falseNegative = 0;

    int n = 0;

    // READ FILE

    // get start time
    clock_gettime(CLOCK_MONOTONIC, &start);

    for (int i = 0; i < 3; i++)
    {
        uniqueWordListLength[i] = readFile(fileNames[i], fileLengths[i], &ppUniqueWordList[i]);
        n += uniqueWordListLength[i];
    }

    // get end time and calculate time taken
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;

    printf("\nTotal unique words from read files: %d\n", n);
    printf("Process time to read file: %lf\n", time_taken);

    // get start time
    clock_gettime(CLOCK_MONOTONIC, &start);

    // m = bloom filter size
    int m = -(n * log(MAX_FP) / pow(log(2), 2));

    // initialise bloom filter
    bool *bloomfilter = (bool *)malloc(m * sizeof(bool));
    for (int i = 0; i < m; i++)
        bloomfilter[i] = false;

    // get end time and calculate time taken
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;

    printf("\nCalculated bit array size for bloom filter: %d\n", m);
    printf("Process time to create bit array of bloom filter: %lf\n", time_taken);

    // INSERT SECTION

    // get start time
    clock_gettime(CLOCK_MONOTONIC, &start);

    // insert into bloom filter
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < uniqueWordListLength[i]; j++)
        {
            insert(bloomfilter, m, ppUniqueWordList[i][j], strlen(ppUniqueWordList[i][j]));
        }
    }

    // get end time and calculate time taken
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;

    printf("\nProcess time to insert all words into bloom filter: %lf\n", time_taken);

    // READ QUERY FILE

    // get start time
    clock_gettime(CLOCK_MONOTONIC, &start);

    readQueryFile(queryName, queryLength, &ppQueryWordList, &pQueryWordTruthList);

    // get end time and calculate time taken
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;

    printf("\nTotal words from reading query file: %d\n", queryLength);
    printf("Process time to read all words from query file: %lf\n", time_taken);

    // QUERY SECTION

    // get start time
    clock_gettime(CLOCK_MONOTONIC, &start);

    // getting results from query
    for (int i = 0; i < queryLength; i++)
    {
        results[i] = query(bloomfilter, m, ppQueryWordList[i], strlen(ppQueryWordList[i]));
    }

    // get end time and calculate time taken
    clock_gettime(CLOCK_MONOTONIC, &end);
    time_taken = (end.tv_sec - start.tv_sec) * 1e9;
    time_taken = (time_taken + (end.tv_nsec - start.tv_nsec)) * 1e-9;

    printf("\nProcess time to query all words: %lf\n", time_taken);

    // calculating results truth table
    for (int i = 0; i < queryLength; i++)
    {
        if (results[i] == 1 && pQueryWordTruthList[i] == 1)
            truePositive++;
        else if (results[i] == 1 && pQueryWordTruthList[i] == 0)
            falsePositive++;
        else if (results[i] == 0 && pQueryWordTruthList[i] == 1)
            falseNegative++;
        else
            trueNegative++;
    }
    printf("\nTotal number of true positive: %d \n", truePositive);
    printf("Total number of false positive: %d \n", falsePositive);
    printf("Total number of false negative: %d \n", falseNegative);
    printf("Total number of true negative: %d \n", trueNegative);

    // clean up
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < uniqueWordListLength[i]; j++)
        {
            free(ppUniqueWordList[i][j]);
        }
        free(ppUniqueWordList[i]);
    }

    // get end time for entire code and calculate time taken
    clock_gettime(CLOCK_MONOTONIC, &endComp);
    time_taken = (endComp.tv_sec - startComp.tv_sec) * 1e9;
    time_taken = (time_taken + (endComp.tv_nsec - startComp.tv_nsec)) * 1e-9;
    printf("\nProcess time of entire code: %lf\n", time_taken);

    return 0;
}

/*
An AP Hash function
adapted from https://www.programmingalgorithms.com/algorithm/ap-hash/c/
*/
uint32_t APHash(char *word, int wordLength)
{
    uint32_t hash = 0xAAAAAAAA;
    for (int i = 0; i < wordLength; word++, i++)
    {
        hash ^= ((i & 1) == 0) ? ((hash << 7) ^ (*word) * (hash >> 3)) : (~(hash << 11) + ((*word) ^ (hash >> 5)));
    }
    return hash;
}

/*
A Jenkins Hash function
adapted from http://www.burtleburtle.net/bob/hash/doobs.html
*/
uint32_t JenkinsHash(char *word, int wordLength)
{
    int i = 0;
    uint32_t hash = 0;
    while (i != wordLength)
    {
        hash += *word;
        hash += hash << 10;
        hash ^= hash >> 6;
        i++;
    }
    hash += hash << 3;
    hash ^= hash >> 11;
    hash += hash << 15;
    return hash;
}

/*
function to read a file and store all unique word into a unique word list
modified from fit3143 COUNT UNIQUE WORDS assignment video brief
*/
int readFile(const char *fileName, int fileLength, char ***pppUniqueWordList)
{
    // read file
    FILE *file = fopen(fileName, "r");

    // initialise variables
    char word[MAX_WORD_LENGTH];
    char **wordList = NULL;
    char **uniqueWordList = NULL;
    int uniqueWordListLength = 0;

    // initialise word list and store all word from file into word list
    wordList = (char **)malloc(fileLength * sizeof(char *));
    for (int i = 0; i < fileLength; i++)
    {
        fscanf(file, "%s", word);
        wordList[i] = strdup(word);
    }

    fclose(file);

    // initialise unique word list
    uniqueWordList = (char **)malloc(fileLength * sizeof(char *));

    // store all unique word from word list into unique word list
    for (int i = 0; i < fileLength; i++)
    {
        // convert to lower case (case-insensitive)
        for (int j = 0; wordList[i][j]; j++)
        {
            wordList[i][j] = tolower(wordList[i][j]);
        }

        // comparison to check for unique words
        if (!isWordInList(wordList[i], uniqueWordList, uniqueWordListLength))
        {
            uniqueWordList[uniqueWordListLength] = strdup(wordList[i]);
            uniqueWordListLength++;
        }
    }

    // reallocating memory for unique word list
    uniqueWordList = realloc(uniqueWordList, uniqueWordListLength * sizeof(char *));

    *pppUniqueWordList = uniqueWordList;
    return uniqueWordListLength;
}

/*
function to read a file and store all unique word into a unique word list
modified from fit3143 COUNT UNIQUE WORDS assignment video brief
*/
void readQueryFile(const char *fileName, int fileLength, char ***pppQueryList, int **ppQueryTruthList)
{
    // read file
    FILE *file = fopen(fileName, "r");

    // initialise variables
    char word[MAX_WORD_LENGTH];
    int exist;
    char **queryList = NULL;
    int *queryTruthList = NULL;
    // int queryListLength = 0;

    // initialise word list and store all word from file into word list
    queryList = (char **)malloc(fileLength * sizeof(char *));
    queryTruthList = (int *)malloc(fileLength * sizeof(int));
    for (int i = 0; i < fileLength; i++)
    {
        fscanf(file, "%s %d", word, &exist);
        queryList[i] = strdup(word);
        queryTruthList[i] = exist;
    }

    fclose(file);

    *pppQueryList = queryList;
    *ppQueryTruthList = queryTruthList;
}

/* function to check if word exist in word list*/
bool isWordInList(char *word, char **wordList, int wordListLength)
{
    for (int i = 0; i < wordListLength; i++)
    {
        if (strcmp(word, wordList[i]) == 0)
            return true;
    }
    return false;
}

/* function to insert word into bloomfilter after hashing with simpleAPHash */
void insert(bool *bloomfilter, int filterSize, char *word, int wordLength)
{
    uint32_t hash1 = APHash(word, wordLength);
    uint32_t hash2 = JenkinsHash(word, wordLength);
    bloomfilter[hash1 % filterSize] = true;
    bloomfilter[hash2 % filterSize] = true;
}

int query(bool *bloomfilter, int filterSize, char *word, int wordLength)
{

    uint32_t hash1 = APHash(word, wordLength);
    uint32_t hash2 = JenkinsHash(word, wordLength);
    // if (bloomfilter[hash1 % filterSize])
    // {
    //     return 1;
    // }
    if (bloomfilter[hash1 % filterSize] && bloomfilter[hash2 % filterSize])
    {
        return true;
    }
    return 0;
}
