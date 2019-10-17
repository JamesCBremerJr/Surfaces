/* 
 * dstrip.c
 *
 * A small utility to strip the debugging header from fortran files.  In fact,
 * dstrip reads from a file specified on the command line, reads it until
 * a line with 50 or more "c"s appears, and then copies the remainder of the
 * file to standard out.
 */
#include <stdio.h>

#define N 10000

int countcc(char *str)
/*
 * Return the number of "c" characters in a string.
 */
{
	int count = 0;
	while(*str) {
		if(*str=='c' || *str=='C') 
			count++;
		str++;		
	}
	return count;
}


int main(int argc, char **argv)
{
	FILE *ip;		
	int over_header = 0;
	char line[N+1];

	/* process command line arguments */
	if(argc!=2) {
		printf("dstrip: invalid usage; the correct syntax is: ");
		printf("dstrip <fortran source>\n");
		return -1;
	}

	/* open the input file */
	ip = fopen(argv[1], "r");
	if(ip==NULL) {
		printf("dstrip: unable to open input file %s\n", argv[1]);
		return;
	}

	/* read the input file line by line */
	while(fgets(line, N, ip)!=NULL) {
		if(countcc(line) > 50) 
			over_header=1;
		
		if(over_header)
		        fprintf(stdout, "%s",line);
		else
			fprintf(stdout, "\n");
	}
	
	/* if we have gotten to the end of the file without hitting the debug header,
      then just output the whole file */
	if(!over_header) {
		fseek(ip, 0, SEEK_SET);
		while(fgets(line, N, ip)!=NULL)
		  fprintf(stdout, "%s",line);
	}

	/* close the input file */
	fclose(ip);
	return 0;
}
