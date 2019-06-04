
void implicit_convertion_test() {
	// implicit convertion done! Only allowed in C
	int * test = malloc(3 * sizeof(int));
	test[0] = 5;
	printf("%d\n", test[0]);
	free(test);
}