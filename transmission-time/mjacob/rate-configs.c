#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "rate-configs.h"
#include "rates.h"

static const int dbg = 0;

#if 0
  for (k=0; k<WIDTH_COUNT; k++) {
    for (j=0; j<GUARD_COUNT; j++) {
      for (i=0; i<MCS_COUNT; i++) {

        sprintf(filename, "%s/%d-%d-%d.dat", dir, i, j, k);
#endif

#if 0
struct rate_config_info {
  char filename[MAX_RATES];
  char rate_str[];
  int index;
};
#endif

int
ratestr_to_filename(char *ratestr, char *filename) {
  char mcs[10];
  char guard[10];
  char width[10];
  int strlen = 0;

  char *tmp = ratestr;

  strlen = 5;
  strncpy(mcs, tmp, strlen);
  mcs[strlen] = '\0';
  tmp += (strlen+1);

  strlen = 2;
  strncpy(guard, tmp, strlen);
  guard[strlen] = '\0';
  tmp += (strlen+1);

  strlen = 3;
  strncpy(width, tmp, strlen);
  width[strlen] = '\0';
  tmp += (strlen+1);

  if (dbg) printf("mcs = [%s] guard = [%s] width = [%s]\n", mcs, guard, width);

  sprintf(filename, "%d-%d-%d.dat", mcs_str_to_index(mcs), 
      guard_str_to_index(guard),
      width_str_to_index(width));

  return 1;
}

// returns a pointer to an array of strings
char **read_rate_configs(char *filename, int *count)
{
   int i = 0;
   int rc = 0;
   FILE *fp;
   char str[80];
   int rate_count = 0;
   fp = fopen(filename, "r");
   assert(fp);
   
   char **rate_configs = malloc(sizeof(char *) * MAX_RATES);

   for (i=0; i<MAX_RATES; i++) {
     rate_configs[i] = NULL;
   }

   while ((rc = fscanf(fp, "%s", str)) == 1) {
     if (dbg) printf("input = %s\n", str);
     rate_configs[rate_count] = malloc(strlen(str) + 1);
     strcpy(rate_configs[rate_count],str);
     rate_count++;
   }
   if (dbg) printf("rc = %d input = %s\n", rc, str);

   *count = rate_count;

   return rate_configs;
}

// for (k=0; k<WIDTH_COUNT; k++) {
// for (j=0; j<GUARD_COUNT; j++) {
// for (i=0; i<MCS_COUNT; i++) {
// sprintf(filename, "%s/%d-%d-%d.dat", dir, i, j, k);

int file_config_name_to_index(char *filename)
{
  char tmp[80];
  int rc = 0;
  int mcs = 0;
  int guard = 0;
  int width = 0;
  int universal_index = 0;

  rc = sscanf(filename, "%d-%d-%d.dat", &mcs, &guard, &width);
  assert(rc == 3);
  if (rc != 3) {
    printf("Error: Filename not using expected format filename = %s rc = %d\n", filename, rc);
    exit(1);
  }

  strcpy(tmp, mcs_to_str(mcs));
  strcat(tmp, "-");
  strcat(tmp, guard_to_str(guard));
  strcat(tmp, "-");
  strcat(tmp, width_to_str(width));
  strcat(tmp, "=??");

#if 0
  if (universal_rates[1] == NULL) {
    universal_rates_init();
  }
#endif

  if (dbg) printf("mcs = %d guard = %d width = %d tmp = %s\n",
                   mcs, guard, width, tmp);

  universal_index = ratestr_to_universal_index(tmp);
  return universal_index;
}

bool filename_in_rate_config(char *filename, char **rate_configs)
{
  int i = 0;
  for (i=0; i<MAX_RATES; i++) {
    if (rate_configs[i]) {
      if (strcmp(rate_configs[i], filename) == 0) {
        return true;
      }
    }
  }
  return false;
}


void 
print_rate_info(int index)
{
  char ratestr[80];
  char filename[80];
  assert(index >= 1);
  assert(index <= RATES);
  rateindex_to_str(index, ratestr);
  ratestr_to_filename(ratestr, filename);
  printf("Rate_index = %4d %20s %12s", index, ratestr, filename);
}
