void diff_c(double* cumPixels, double* diff, int n){
    int i;

    diff[0] = cumPixels[0];
    for(i = 1; i < n; i++){
        diff[i] = cumPixels[i] - diff[i-1];
    }
    return ;
}

void seg_c(double* diff, double* boxing, int n){
	int i;

	for(i = 1 ; i < n-1; i++){
		boxing[i] = 1;
		if( ((diff[i] < diff[i-1]) && (diff[i] < diff[i+1])) || ((diff[i] > diff[i-1]) && (diff[i] > diff[i+1]))){
			boxing[i] = 0;
		}
	}
	return ;
}

void eros_c(double* boxing, double* blocking, int n){
	int i;
	for(i= 0; i < n-1;i++){
		blocking[i] = boxing[i];
		if ( boxing[i-1] == 0 && boxing[i] == 1 && boxing[i+1] == 0 ){
			blocking[i] = 0;
		}
	}

	for(i = 0; i < n; i++){
		boxing[i] = blocking[i];
	}
	for(i = 1; i < n-1; i++){
		if(blocking[i-1] == 0 && blocking[i]== 1){
			boxing[i-1] = 1;
		}
		if(blocking[i]==1 && blocking[i+1]==0){
			boxing[i+1] = 1;
		}
	}
	return;
}