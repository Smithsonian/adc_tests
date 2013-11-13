#include <stdio.h>
#include <sys/file.h>
#include <string.h>
#include <math.h>
#define USE_MAIN 1

/* The ogp array contains o, g and p of core A followed by the same for cores
 * B, C and D.
 * In this case, since phase can not be determined, each p will be 0.
 * Overlaod_cnt will contain a count of -128 and 127 codes for each core.
 */
typedef struct {
  float ogp[12];
  float avz;
  float avamp;
  int overload_cnt[4];
} og_rtn;

/* Cores are read out in the sequence A, C, B, D , but should be in the natural
 * order in the ogp array. */
int startpos[] = {0, 2, 1, 3};
og_rtn og_from_noise(int len, char *snap) {
  og_rtn rtn;
  float avg, amp;
  int sum, i, cnt, core, code, start_i;

  memset((char *)&rtn, 0, sizeof(rtn));
  for(core = 0; core < 4; core++) {
    start_i = startpos[core];
    cnt = 0;
    sum = 0;
    amp = 0;
    for(i=start_i; i < len; i+= 4) {
      cnt++;
      code = snap[i];
      sum += code;
      if(code == -128 || code == 127) {
        rtn.overload_cnt[core]++;
      }
    }
    avg = (float)sum / cnt;
    for(i=start_i; i < len; i+= 4) {
      amp += fabs((float)snap[i] - avg);
    }
    avg *= (-500.0/256);
    rtn.avz += avg;
    rtn.avamp += amp/cnt;
    rtn.ogp[core*3] = avg;
    rtn.ogp[core*3 + 1] = amp/cnt;
  }
  rtn.avz /= 4;
  rtn.avamp /= 4;
  for(core = 0; core < 4; core++) {
    i = core*3 + 1;
    rtn.ogp[i] = 100 * (rtn.avamp - rtn.ogp[i])/rtn.avamp;
  }
  return(rtn);
}

#if USE_MAIN
int main(int argc, char *argv[]) {
  og_rtn rtn;
  char buf[16384];
  int len;
  FILE * fd;
  char fname[] = "t.og_noise";

  if((fd = fopen(fname, "r")) == NULL) {
    printf("Can't open %s\n", fname);
    return(1);
  }
  memset(buf, 0, sizeof(buf));
  for(len = 0; len < (int)sizeof(buf); len++) {
    if(fscanf(fd, "%hhd", buf+len) == EOF) {
      printf("at %d, found EOF\n", len);
      break;
    }
  }
  printf("%d  %d %d %d %d\n", len, buf[0], buf[1], buf[2], buf[3]);
  rtn = og_from_noise(len, buf); 
  printf("      Offset(mV) Gain(%%) Overflows (adj by .4, .14)\n");
  printf("All    %7.4f  %7.4f\n", rtn.avz, rtn.avamp);
  printf("Core A %7.4f  %7.4f  %3d\n", rtn.ogp[0], rtn.ogp[1], rtn.overload_cnt[0]);
  printf("Core B %7.4f  %7.4f  %3d\n", rtn.ogp[3], rtn.ogp[4], rtn.overload_cnt[1]);
  printf("Core C %7.4f  %7.4f  %3d\n", rtn.ogp[6], rtn.ogp[7], rtn.overload_cnt[2]);
  printf("Core D %7.4f  %7.4f  %3d\n", rtn.ogp[9], rtn.ogp[10], rtn.overload_cnt[3]);
}
#endif
