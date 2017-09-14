#include <stdio.h>
#include <stdlib.h>
#include <curses.h>
#include <signal.h>
#include <string.h>

/* variables globales */

SCREEN *ecran;
FILE *ttyio;
int insert=0;
int xposi=22;
int scroll_x,scroll_y;

/**********************************************************************
  sig_catch : Routine servant a attraper le caractere CTRL-C.
              ou le caractere d'interuption.
 **********************************************************************/

void sig_catch(int scrap)
{
	endwin();
	printf(" CTRL-C detecte \n");
	exit(1);
}


/**********************************************************************

   vmenu menu pour M. Valin

   Ecrit par W. Richards CIDU.
   date : 91/03/20

   fonction : a partir des variables entrees titre, cle, val et help,
	      afficher un ecran sur "/dev/tty" qui permettra l'entree
	      interactive de donnee a meme un menu.

 **********************************************************************/

int vmenu(titre,cle,val,help,nbcle,nfois)
char *titre;
char **cle,**val,**help;
int nbcle,nfois;
{
 int center;
 int i=0;
 char nom_terminal[50];
 char *newlpt;

#if defined (HP)
 initscr();
#endif

 if((ttyio = fopen("/dev/tty","r+")) == (FILE *) NULL)
     return(-1);

 strcpy(nom_terminal,getenv("TERM"));
 ecran = newterm(nom_terminal,ttyio,ttyio);
  refresh();

 (void) signal(SIGINT,sig_catch);
 cbreak();
/* nonl();*/
 keypad(stdscr,TRUE);
 idlok(stdscr,FALSE);
 scrollok(stdscr,FALSE);
 noecho();

 /* Calcul de la variable center  et initialisation de scroll_x et scroll_y */
 center = (COLS - strlen(titre))/2;
 scroll_x=0; scroll_y=0;
 refresh();
 attrset(A_REVERSE);
 mvaddstr(0,center,titre);
 attrset(A_BOLD);
 mvaddstr(LINES-5,0,"--------------------------------------");
 if (insert)
 mvaddstr(LINES-5,COLS-10,"Insert   ");
 else
 mvaddstr(LINES-5,COLS-10,"Overwrite");
 if(!nfois)
 {
      newlpt = strchr(help[0],'\\');
      if((newlpt != (char *) NULL) && (*(newlpt+1) == 'n'))
      {
            *newlpt = ' ';
            *(newlpt+1) = '\n';
      }
      mvaddstr(LINES-4,0,help[0]);
 refresh();
 }
 else
 {
      mvaddstr(LINES-4,0,"s.v.p. donnez une valeur aux clefs non initialisees");
 refresh();
 }
 
 attrset(A_REVERSE);
 mvaddstr(LINES-2,0,"^X:End");
 mvaddstr(LINES-1,0,"^C:Abt");
 mvaddstr(LINES-2,8,"^U:Undo");
 mvaddstr(LINES-1,8,"^A:Help");
 mvaddstr(LINES-2,17,"^O:Over/Ins");
 mvaddstr(LINES-1,17,"^E:End-line");
 mvaddstr(LINES-2,30,"^D:Delete ");
 mvaddstr(LINES-1,30,"^L:Refresh");
 refresh();
 if(nbcle > (LINES - 6))
 {
      mvaddstr(LINES-2,42,"^F:PGdn");
      mvaddstr(LINES-1,42,"^B:Pgup");
      refresh();
 }
 attrset(0);
 /* Impression des champs et valeurs de la page initiale */
 imprimer_ecran(nbcle,cle,val);

 for (i=scroll_y+1;i<nbcle;)
   {
   /* get_chaine retourne +1     pour aller a la prochaine ligne,
			  -1     pour aller a la precedente,
			  1000   pour indiquer CTRL-X ( fin),
			  2000   pour indiquer CTRL-F (page down),
			  3000   pour indiquer CTRL-B (page up),
   */
   i = i + get_chaine(cle,val,help,i,nbcle);
   if (i-scroll_y==0) {
      if (scroll_y) {  /* Scroll up */
	 scroll_y--;
	 move(LINES-6,0);deleteln();
	 move(1,0);insertln();
	 attrset(A_BOLD);
         mvaddstr(i-scroll_y,2,cle[i]);
         mvaddstr(i-scroll_y,20,":\0"); 
         attrset(0);
         mvaddstr(i-scroll_y,22,&val[i][scroll_x]);
         refresh();
	 }
      else {   /* Va a la derniere position */
         i=nbcle-1;
	 scroll_y = nbcle - (LINES - 5);
	 if (scroll_y < 0) scroll_y=0;
	 imprimer_ecran(nbcle,cle,val);
	 }
   }
   if (i==nbcle) { /* Va a la premiere position */
      scroll_y=0;
      imprimer_ecran(nbcle,cle,val);
      i=1;
      };
   if (i-scroll_y==LINES-5) {  /* Scroll down */
      scroll_y++;
      move(1,0);deleteln();
      move(LINES-6,0);insertln();
      attrset(A_BOLD);
      mvaddstr(i-scroll_y,2,cle[i]);
      mvaddstr(i-scroll_y,20,":\0"); 
      attrset(0);
      mvaddstr(i-scroll_y,22,&val[i][scroll_x]);
      refresh();
      }
    if (i > 3000) { /* CTRL-B page up */
       if (nbcle > (LINES - 6) ) {
          if (scroll_y - LINES + 7 > 0) 
	      scroll_y = scroll_y - LINES + 7;
          else scroll_y = 0;
          if  (i < 3000 + scroll_y + LINES - 7 && i > 3000 + scroll_y ) ;
          else i = scroll_y + 1 + 3000;

       imprimer_ecran(ecran,cle,val);
       }
       i = i - 3000;
       xposi = 22; scroll_x = 0;
    }
    if (i > 2000) { /* CTRL-F page down */
       if (nbcle > (LINES - 6) ) {
          if (2*(LINES - 7) + scroll_y <= nbcle) 
	      scroll_y = scroll_y + LINES - 7;
          else scroll_y = nbcle - LINES + 5;
          if  (i < 2000 + scroll_y + LINES - 7 && i > 2000 + scroll_y ) ;
          else i = scroll_y + 1 + 2000;

       imprimer_ecran(ecran,cle,val);
       }
       i = i - 2000;
       xposi = 22; scroll_x = 0;
    }
   }
 endwin(); 
 fclose(ttyio);
 return(1);
}


/**********************************************************************
  get_chaine : Routine principale servant a obtenir la chaine de
	       caractere finale. Cette routine sert aussi a gerer
	       les caracteres speciaux d'entree-sortie.
 **********************************************************************/
int get_chaine(cle,chaine,help,num,nbcle)
char **cle;
char **chaine;
char **help;
int num,nbcle;
{
char tmp[500];
char blank[130];
char *ptr, *newlpt;
int key,i,j,x,y,itmp;

for (i=0;i<COLS-22&&i<130;i++) blank[i]=' ';
for (i=0;i<500;i++) tmp[i]='\0';
move(num-scroll_y,xposi);
refresh();
strcpy(tmp,chaine[num]);
ptr = &tmp[xposi-22+scroll_x];

for (i=xposi-22+scroll_x; ;i++)
{
  key=getch();
  switch (key)
    {

       case 4 : /* CTRL-D effacer un caractere sous le curseur  */
	       getyx(stdscr,y,x);
	       if (strlen(tmp) > scroll_x+x-22) 
               {
	          if (strlen(tmp) > scroll_x+COLS-2-22)
	              mvinsch(y,COLS-2,tmp[scroll_x+COLS-22-2]);
                  move(y,x);
	          delch();
                  refresh();
	          for (j=0;j<strlen(tmp)+1;j++) ptr[j]=ptr[j+1];
               }
	       break;
       case 5 :  /* si CTRL-E on deplace le curseur a la fin du champ */
	       getyx(stdscr,y,x);
	       itmp = strlen(tmp);
	       if (itmp > scroll_x && itmp < COLS-2-22+scroll_x){
		    move(y,x+itmp-scroll_x -x +22);
                   } 
               else  if (itmp < COLS -2 -22)
		         {scroll_x = 0;
			  imprimer_ecran(nbcle,cle,chaine);
			  move(y,22+itmp);
                         } 
                      else {scroll_x = itmp-COLS+2+22;
			    imprimer_ecran(nbcle,cle,chaine);
			    mvaddstr(y,22,&tmp[itmp-COLS+2+22]);
                         } 
               refresh();
	       i = itmp;
	       ptr = &tmp[itmp];
	       break;
       case KEY_BACKSPACE :/*si Backspace, on fait backspace mais on          */
       case 8 :            /*peu aussi effacer le caractere prec. si overwrite*/
	       if (insert) {
	           i = i -2;
	           if (i<-1) i++;
	           else {
			ptr--;
			itmp=strlen(ptr);
			for (j=0;j<itmp+1;j++) ptr[j]=ptr[j+1];
			getyx(stdscr,y,x);
			if (strlen(tmp) > scroll_x+COLS-2-22)
			   { move(y,COLS-1);
			   printw("%c",tmp[scroll_x+COLS-22-2]);}
			if (x == 22) {
			    move_right(nbcle,cle,chaine); delch();}
			else mvdelch(y,x-1);
			}
                    refresh();
                   break;}
		      
       case KEY_LEFT : /* Si key_left, on deplace le curseur a gauche */
		       /* ATTENTION , si Backspace et pas en insert mode
			  on passe a KEY_LEFT */
	       i = i -2;
	       getyx(stdscr,y,x);
	       if ((x == 22) && scroll_x) {
		  move_right(nbcle,cle,chaine);
		  if (strlen(chaine[num])<=scroll_x)
		     insch(tmp[scroll_x]);
		  else printw("%c",tmp[scroll_x]);
		  move(y,x);
		  ptr--;}
	       else if (i<-1) i++;
	            else {printw("\b");ptr--;}
               refresh();
               break;

       case KEY_UP : /* Si key_up, on sauvegarde la chaine temporaire et on remonte */
	       compare(tmp,chaine,num);
	       getyx(stdscr,y,xposi);
	       return(-1);

       case KEY_RIGHT : /* Si key_right, on deplace le curseur a droite */
	       getyx(stdscr,y,x);
	       ptr++;
	       if (x == COLS -2) {
		  move_left(nbcle,cle,chaine);
                  printw("%c",ptr[0]);
                  move(y,x);}
	       else move(y,x+1);
               refresh();
	       break;

       case KEY_DOWN : /* Si key_down, on sauvegarde la chaine temporaire et on descend */
       case 9 :
	       compare(tmp,chaine,num);
	       getyx(stdscr,y,xposi);
	       return(1);

       case 343 : /* Si Carriage_return, on sauvegarde la chaine temp. et on descend */
       case 10 :
	       compare(tmp,chaine,num);
	       xposi=22;
	       if (scroll_x) {
	           scroll_x=0;
	           imprimer_ecran(nbcle,cle,chaine);}
	       return(1);

       case 21 : /* Si CTRL-U, on fait un UNDO pour la ligne courrante */
	       if (scroll_x) {
	           scroll_x=0;
	           imprimer_ecran(nbcle,cle,chaine);}
	       strcpy(tmp,chaine[num]);
	       i= -1;
	       ptr=tmp;
	       mvaddstr(num,22,blank);
	       mvaddstr(num,22,tmp);
	       move(num,22);
               refresh();
	       break;

       case 24 : /* Si CTRL-X, on fait une SORTIE normal */
	       compare(tmp,chaine,num);
	       return(1000);

       case 6  : /* Si CTRL-F, on veut faire un page down */
	       compare(tmp,chaine,num);
	       return(2000);

       case 2  : /* Si CTRL-B, on veut faire un page up   */
	       compare(tmp,chaine,num);
	       return(3000);

       case 1 : /* Si CTRL-A, on imprime le Help, (ASSIST) pour le champ courrant */
	       getyx(stdscr,y,x);
               for (j=LINES-4;j<LINES-2;j++)
               {
               mvaddstr(j,0,blank);
               mvaddstr(j,22,blank);
               }
	       mvaddstr(LINES-5,0,blank);
	       attrset(A_REVERSE);
	       mvaddstr(LINES-5,0,"HELP : ");
	       mvaddstr(LINES-5,7,cle[num]);
	       attrset(A_BOLD);
               newlpt = strchr(help[num],'\\');
               if((newlpt != (char *) NULL) && (*(newlpt+1) == 'n'))
               {
                       *newlpt = ' ';
                       *(newlpt+1) = '\n';
               }
               move(LINES-4,0);
               printw("%s",help[num]);
	       attrset(0);
	       move(y,x);
	       i--;
               refresh();
               break;

       case 15 : /* Si CTRL-O, on toggle entre Overwrite et Insert mode */
	       getyx(stdscr,y,x);
	       attrset(A_BOLD);
               if (insert) 
	         mvaddstr(LINES-5,COLS-10,"Overwrite");
	       else
	         mvaddstr(LINES-5,COLS-10,"Insert   ");
	       attrset(0);
	       insert = ! insert;
	       move(y,x);
	       i--;
               refresh();
	       break;

       case 12  : /* Si CTRL-L, on Retrace l'ecran */
/*		  clearok(stdscr,TRUE); */
		  refresh();
/*		  clearok(stdscr,FALSE); */
		  i--;
		  break;

       default : /* Par default, on met a jour la chaine */
	       if ((key<32) || (key>126)) i--;
               refresh();
	       break;
     }
if ((key >=32) && (key<=126)) /* Mise a jour de la chaine pour caractere imprimable */
{                             /* et normaux */
  getyx(stdscr,y,x);
  /* au cas ou on a deplace le curseur avec les fleches et on veut inserer */
  if (x-22+scroll_x > strlen(tmp))
      { for (j=strlen(tmp)-1;j<x-22+scroll_x; j++) tmp[j]=' '; tmp[j]='\0';};
  if (x == COLS -2) {
     move_left(nbcle,cle,chaine);
     move(y,x-1);
     printw("%c",key);
     printw("%c",ptr[0]);
     move(y,x);
     }
  else if (insert) {
           insch(key);
           getyx(stdscr,y,x);
           mvaddstr(y,COLS-1," ");
           move(y,x+1);
           }
         else 
           printw("%c",key);
  if (insert) {
     itmp=strlen(ptr);
     for (j=itmp;j>-1;j--) ptr[j+1]=ptr[j];
     *ptr++ = (char) key;}
  else *ptr++ = (char) key;

}
refresh();
}
}


/**********************************************************************
  compare : Routine servant a comparer l'entree d'un usager avec
	    l'entree de default donne par la liste de depart.

 **********************************************************************/
compare(tmp,chaine,num)
char *tmp,**chaine;
int num;
{
 char *tmpp;
 if (strcmp(tmp,chaine[num])) 
   {
     free(chaine[num]);
     chaine[num] = (char *) malloc(strlen(tmp)+1);
     strcpy(chaine[num],tmp);
   }
}


imprimer_ecran(nbcle,cle,val)
int nbcle;
char **cle,**val;

{
int i;
char tmp[30];

for ( i=scroll_y+1;((i<nbcle) && (i-scroll_y<LINES-5)); i++) 
  {
   attrset(A_BOLD);
   move(i-scroll_y,0);deleteln();insertln();
   mvaddstr(i-scroll_y,2,cle[i]);
   mvaddstr(i-scroll_y,20,":\0"); 
   attrset(0);
   sprintf(tmp,"%%.%-ds",COLS-22-2);
   mvprintw(i-scroll_y,22,tmp,&val[i][scroll_x]);
   mvaddstr(i-scroll_y,COLS-1," ");
  }
/* clearok(stdscr,TRUE); */
 refresh();
/* clearok(stdscr,FALSE); */
}


move_left(nbcle,cle,val)
int nbcle;
char **cle,**val;

{
int i,x,y;

getyx(stdscr,y,x);
for ( i=scroll_y+1;((i<nbcle) && (i-scroll_y<LINES-5)); i++) 
  {
   mvdelch(i-scroll_y,22);
   if (strlen(val[i]) > COLS-2-22+scroll_x) {
       move(i-scroll_y,COLS-2);
       printw("%c",val[i][scroll_x+COLS-1-22]);};
  }
move(y,x);
scroll_x++;
}


move_right(nbcle,cle,val)
int nbcle;
char **cle,**val;

{
int i,x,y;

if (scroll_x) {
   scroll_x--;
   getyx(stdscr,y,x);
   for ( i=scroll_y+1;((i<nbcle) && (i-scroll_y<LINES-5)); i++) 
     {
          move(i-scroll_y,22);
	  if (strlen(val[i]) > scroll_x+1)
              insch(val[i][scroll_x]);
          mvaddstr(i-scroll_y,COLS-1," ");
     }
   move(y,x);
}
}
