#if !defined(F2Cl)
#define F2Cl int
#endif
int fnom_index(int iun);
int error_msg(char *function_name, int errcode, int errlevel);
int file_index(int iun);
void build_burp_prim_keys(burp_record *brpk, word *keys,
                                 burp_record *mask, word *mskkeys,
				 int index, int mode);
void build_burp_info_keys(word *buf, word *keys, int index, int mode);
void build_fstd_info_keys(word *buf, word *keys, int index, int mode);
void build_fstd_prim_keys(word *buf, word *keys, word *mask, word *mskkeys,
                                int index, int mode);

int c_qdfmsig(int iun, char* newappl);
ftnword f77name(qdfmsig)(ftnword *fiun,char *appl,F2Cl l1);
ftnword f77name(mrbdel)(word *buf, ftnword *f_number);
int c_mrbdel(void *buffer, int number);
ftnword f77name(fstvoi)(ftnword *f_iun,char *options,F2Cl l1);
int c_fstvoi(int iun,char *options);
ftnword f77name(fstouv)(ftnword *f_iun, char *options, F2Cl l1);
int c_fstouv(int iun, char *options);
ftnword f77name(secateur)(char *filename, ftnword *f_where, F2Cl l1);
void c_fst_env_var(char *cles, int index, char *content);
