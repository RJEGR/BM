Lets figure out how to make index for ASVS/otus sequence in order to define as sequence id (based on the callahan proposal to call ASVs by non-character-numeric identifier):

```c
# include <stdio.h>
  int main(){
   printf("Hello World!");
   return 0;
}
```

save in a file of name `helloWorld.c` and compile with:

1) `clang helloWorld.c -o helloWorld`

2) `gcc helloWorld.c`

Lets index a sequence from ncbi with header fasta of `>LT840348.1 Knoutsodonta sp. B GF-2017 mitochondrial partial COI gene for cytochrome oxidase subunit 1, specimen voucher BAU02983`:

```c
// C program to find Burrows Wheeler transform of 
// a given text 
#include<stdio.h> 
#include<string.h> 
#include<stdlib.h> 

// Structure to store data of a rotation 
struct rotation 
{ 
	int index; 
	char *suffix; 
}; 

// Compares the rotations and 
// sorts the rotations alphabetically 
int cmpfunc (const void *x, const void *y) 
{ 
	struct rotation *rx = (struct rotation *)x; 
	struct rotation *ry = (struct rotation *)y; 
	return strcmp(rx -> suffix, ry -> suffix); 
} 

// Takes text to be transformed and its length as 
// arguments and returns the corresponding suffix array 
int *computeSuffixArray(char *input_text, int len_text) 
{ 
	// Array of structures to store rotations and 
	// their indexes 
	struct rotation suff[len_text]; 

	// Structure is needed to maintain old indexes of 
	// rotations after sorting them 
	for(int i = 0; i < len_text; i++) 
	{ 
		suff[i].index = i; 
		suff[i].suffix = (input_text+i); 
	} 

	// Sorts rotations using comparison function defined above 
	qsort(suff, len_text, sizeof(struct rotation), cmpfunc); 

	// Stores the indexes of sorted rotations 
	int *suffix_arr = (int *) malloc (len_text * sizeof(int)); 
	for (int i = 0; i < len_text; i++) 
		suffix_arr[i] = suff[i].index; 

	// Returns the computed suffix array 
	return suffix_arr; 
} 

// Takes suffix array and its size as arguments and returns 
// the Burrows - Wheeler Transform of given text 
char *findLastChar(char *input_text, int *suffix_arr, int n) 
{ 
	// Iterates over the suffix array to find 
	// the last char of each cyclic rotation 
	char *bwt_arr = (char *) malloc (n * sizeof(char)); 
	int i; 
	for (i = 0; i < n; i++) 
	{ 
		// Computes the last char which is given by 
		// input_text[(suffix_arr[i] + n - 1) % n] 
		int j = suffix_arr[i] - 1; 
		if (j < 0)	 
			j = j + n; 
		
		bwt_arr[i] = input_text[j]; 
	} 

	bwt_arr[i] = '\0'; 

	// Returns the computed Burrows - Wheeler Transform 
	return bwt_arr; 
} 

// Driver program to test functions above 
int main() 
{ 
	char input_text[] = "GGAACTGGGTTGAGATTACTCATTCGTTTTGAATTAGGAACAGCAGGAGCTTTTTTAGGAGATGATCATTTTTATAATGTAATTGTTACGGCTCATGCCTTTGTTATAATTTTTTTTATGGTAATACCTCTAATAATTGGAGGATTTGGTAACTGAATAGTACCATTATTAATTGGTGCACCAGATATAAGTTTTCCTCGAATAAACAATATAAGGTTTTGGTTACTACCTCCATCTTTTATTCTACTCTTATGTTCTACTTTGATAGAAGGAGGAGCAGGAACAGGCTGAACTGTTTACCCTCCTTTATCTGGTCCAGTAGGCCATAGAGGGGCTTCGGTAGATCTGGCTATTTTTTCTTTACATTTGGCGGGTGCATCTTCATTGTTAGGAGCAATCAATTTTATTACAACTATTTTTAATATACGGTCGCCAGCAATAAGTATGGAACGTGTCAGGTTATTTGTTTGATCTGTCTTAGTGACTGCTTTTTTATTATTATTATCTTTACCTGTTTTGGCAGGAGCAATTACTATACTATTAACAGATCGTAATTTTAACACTAGTTTTTTCGACCCGGCAGGAGGAGGAGATCCAATTTTGTACCAAC$"; 
	int len_text = strlen(input_text); 

	// Computes the suffix array of our text 
	int *suffix_arr = computeSuffixArray(input_text, len_text); 

	// Adds to the output array the last char of each rotation 
	char *bwt_arr = findLastChar(input_text, suffix_arr, len_text); 

	printf("Input text : %s\n", input_text); 
	printf("Burrows - Wheeler Transform : %s\n", bwt_arr); 
	return 0; 
}
```

Save as `BWT.c` and compile:

1) `clang BWT.c -o BWT`

2) `./BWT`

Lets return the inverse of a sequence indexed:

```C
// C program to find inverse of Burrows 
// Wheeler transform 
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h> 

// Structure to store info of a node of 
// linked list 
struct node { 
	int data; 
	struct node* next; 
}; 

// Compares the characters of bwt_arr[] 
// and sorts them alphabetically 
int cmpfunc(const void* a, const void* b) 
{ 
	const char* ia = (const char*)a; 
	const char* ib = (const char*)b; 
	return strcmp(ia, ib); 
} 

// Creates the new node 
struct node* getNode(int i) 
{ 
	struct node* nn = 
		(struct node*)malloc(sizeof(struct node)); 
	nn->data = i; 
	nn->next = NULL; 
	return nn; 
} 

// Does insertion at end in the linked list 
void addAtLast(struct node** head, struct node* nn) 
{ 
	if (*head == NULL) { 
		*head = nn; 
		return; 
	} 
	struct node* temp = *head; 
	while (temp->next != NULL) 
		temp = temp->next; 
	temp->next = nn; 
} 

// Computes l_shift[] 
void* computeLShift(struct node** head, int index, 
					int* l_shift) 
{ 
	l_shift[index] = (*head)->data; 
	(*head) = (*head)->next; 
} 

void invert(char bwt_arr[]) 
{ 
	int i,len_bwt = strlen(bwt_arr); 
	char* sorted_bwt = (char*)malloc(len_bwt * sizeof(char)); 
	strcpy(sorted_bwt, bwt_arr); 
	int* l_shift = (int*)malloc(len_bwt * sizeof(int)); 

	// Index at which original string appears 
	// in the sorted rotations list 
	int x = 4; 

	// Sorts the characters of bwt_arr[] alphabetically 
	qsort(sorted_bwt, len_bwt, sizeof(char), cmpfunc); 

	// Array of pointers that act as head nodes 
	// to linked lists created to compute l_shift[] 
	struct node* arr[128] = { NULL }; 

	// Takes each distinct character of bwt_arr[] as head 
	// of a linked list and appends to it the new node 
	// whose data part contains index at which 
	// character occurs in bwt_arr[] 
	for (i = 0; i < len_bwt; i++) { 
		struct node* nn = getNode(i); 
		addAtLast(&arr[bwt_arr[i]], nn); 
	} 

	// Takes each distinct character of sorted_arr[] as head 
	// of a linked list and finds l_shift[] 
	for (i = 0; i < len_bwt; i++) 
		computeLShift(&arr[sorted_bwt[i]], i, l_shift); 

	printf("Burrows - Wheeler Transform: %s\n", bwt_arr); 
	printf("Inverse of Burrows - Wheeler Transform: "); 
	// Decodes the bwt 
	for (i = 0; i < len_bwt; i++) { 
		x = l_shift[x]; 
		printf("%c", bwt_arr[x]); 
	} 
} 

// Driver program to test functions above 
int main() 
{ 
	char bwt_arr[] = "CTCATTGGGCTGGGTTTGCTTGCTCTCGTTTTCCTATAAAAATTCTGTTTTTTATCTTATTAGAATTTCGCTGGCGGGCGTCGTTCGCGACGTCGCATCATATATATTATATTGCAAGAAGGGGTGCTCGCTTATTTTAGATCTTCTAAACCGTAAATCTAACAGAGGTCGACACAGGGGATCCGCTCTTATATAGTGTAAACACTAATGTTTCAGTATATTATTAAAAGGACCCCAGAATATTACTTTGTTCACTGGGGGGTACTTCTGGTGGGGGGGGATATAAATAAGAAATAGAGTCGTGGGGGATAAT$AAAAAAAAATATCATTCAGGACTTTCTCTGTAGAGGTCTAGAATGGTAGGCTGTTGTGTTACTGAATTGAAACGTAATGGATTGGTCATTATTATCCAAGTTTGATCAATTACTTGTTCTTTTTTTCCAGCCTAAGCTCCTGTATCTTAAAGCATAACCTGTTATGGACATTTCCTACTTATGCTTTACTCTATAATTTGGAAATGCTGTATCATAATTGCTTCAAGTTGTTATTTATTAATTTTCGCTTCTTTTTTTTCGTAATTCAATTTTATCGTTGGGAATATTTTCTCGATA"; 
	invert(bwt_arr); 
	return 0; 
} 
```

Save in a file of name `BWTi.c`, and then run:

1) `clang BWTi.c -o BWTi`

2) `./BWTi`

## Reference

1. https://www.geeksforgeeks.org/burrows-wheeler-data-transform-algorithm/
2. https://www.geeksforgeeks.org/inverting-burrows-wheeler-transform/