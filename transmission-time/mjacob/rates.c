#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "rates.h"

// This will be determined at runtime.
int RATES = 0;

static const int dbg = 0;


static char *mcs_strs[] = {
  "1S-I1", "1S-I2", "1S-I3", "1S-I4", "1S-I5", "1S-I6", "1S-I7", "1S-I8",
  "2S-I1", "2S-I2", "2S-I3", "2S-I4", "2S-I5", "2S-I6", "2S-I7", "2S-I8",
  "3S-I1", "3S-I2", "3S-I3", "3S-I4", "3S-I5", "3S-I6", "3S-I7", "3S-I8",
};

static char *guard_strs[] = {
  "LG", "SG"
};

static char *width_strs[] = {
  "20M", "40M",
};

static char *phyrate_strs[] = {
  "6.5", "13", "19.5", "26", "39", "52", "58.5", "65", 
  "13", "26", "39", "52", "78", "104", "117", "130", 
  "19.5", "39", "58.5", "78", "117", "156", "175.5", "195", 
  "7.2", "14.4", "21.7", "28.9", "43.3", "57.8", "65", "72.2",
  "14.4", "28.9", "43.3", "57.8", "86.7", "115.6", "130", "144.4",
  "21.7", "43.3", "65", "86.7", "130", "173.3", "195", "216.7",
  "13.5", "27", "40.5", "54", "81", "108", "121.5", "135",
  "27", "54", "81", "108", "162", "216", "243", "270",
  "40.5", "81", "121.5", "162", "243", "324", "364.5", "405",
  "15", "30", "45", "60", "90", "120", "135", "150",
  "30", "60", "90", "120", "180", "240", "270", "300",
  "45", "90", "135", "180", "270", "360", "405", "450"
};

// This is the list of all 96 rates
static char *universal_rates[MAX_RATES+1] = {NULL};

// This is the list of only the available rates
// If there are 16 rates everything after that will be NULL
static char *available_rates[MAX_RATES+1] = {NULL};

int
mcs_str_to_index(char *mcs) 
{
  int i = 0;
  for (i=0; i<sizeof(mcs_strs) / sizeof(char *); i++) {
    if (strcmp(mcs, mcs_strs[i]) == 0) {
      return i;
    }
  }
  printf("Could not find MCS string\n");
  assert(1 == 0);
  return 0;
}


char *mcs_to_str(int mcs)
{
  assert(mcs >= 0);
  assert(mcs < (sizeof(mcs_strs) / sizeof(char *)));
  return mcs_strs[mcs];
}

int
guard_str_to_index(char *guard) 
{
  int i = 0;
  for (i=0; i<sizeof(guard_strs) / sizeof(char *); i++) {
    if (strcmp(guard, guard_strs[i]) == 0) {
      return i;
    }
  }
  printf("Could not find guard interval string\n");
  assert(1 == 0);
  return 0;
}

char *guard_to_str(int guard)
{
  assert(guard >= 0);
  assert(guard < (sizeof(guard_strs) / sizeof(char *)));
  return guard_strs[guard];
}

int
width_str_to_index(char *width) 
{
  int i = 0;
  for (i=0; i<sizeof(width_strs) / sizeof(char *); i++) {
    if (strcmp(width, width_strs[i]) == 0) {
      return i;
    }
  }
  printf("Could not find channel width string\n");
  assert(1 == 0);
  return 0;
}

char *width_to_str(int width)
{
  assert(width >= 0);
  assert(width < (sizeof(width_strs) / sizeof(char *)));
  return width_strs[width];
}

static void
universal_rateindex_fill_str(int i, char *ratestr)
{
  int index = i-1;
  int mcs = (index % 24);
  int guard = (index / 24) % 2;
  int width = (index / 48);

  if (dbg) printf("index = %d mcs = %d guard = %d width = %d\n", index, mcs, guard, width);

  sprintf(ratestr, "%s-%s-%s=%s", mcs_strs[mcs], guard_strs[guard], width_strs[width], phyrate_strs[index]);

  fflush(stdout);
}

static void
universal_rates_init()
{
  int i = 0; 
  char *ratestr = 0;
  // printf("universal_rates_init\n");
  for (i=1; i<=MAX_RATES; i++) {
     ratestr = malloc(RATE_STR_MAX * sizeof(char));
     assert(ratestr);
     universal_rateindex_fill_str(i, ratestr);
     universal_rates[i] = ratestr;
     if (dbg) printf("universal_rates[%d] = %s\n", i, universal_rates[i]); 
  }
}

bool ratestr_is_in_available_rates(char *ratestr)
{
  int i=0;
  for (i=1; i<=MAX_RATES; i++) {
    if (available_rates[i] != NULL) {
      if (strcmp(ratestr, available_rates[i]) == 0)  {
        return true;
      }
    }
  }
  return false;
}

int ratestr_to_available_index(char *ratestr)
{
  int i=0;
  for (i=1; i<=MAX_RATES; i++) {
    if (available_rates[i] != NULL) {
      if (strcmp(ratestr, available_rates[i]) == 0)  {
        return i;
      }
    }
  }
  assert(1 == 0);
  return false;
}

char *
available_rate_str(int rate) 
{
   assert(available_rates[rate] != 0);
   return(available_rates[rate]);
}

void
available_rate_init(int rate, int mcs, int guard, int width) 
{
  char tmp[RATE_STR_MAX];
  int universal_index = 0;
  tmp[0] = '\0';
   
  strcpy(tmp, mcs_strs[mcs]);
  strcat(tmp, "-");
  strcat(tmp, guard_strs[guard]);
  strcat(tmp, "-");
  strcat(tmp, width_strs[width]);
  strcat(tmp, "=??");

  if (universal_rates[1] == NULL) {
    universal_rates_init();
  }

  if (dbg) printf("rate = %d mcs = %d guard = %d width = %d tmp = %s\n",
                   rate, mcs, guard, width, tmp);

  universal_index = ratestr_to_universal_index(tmp); 
  if (dbg) printf("universal_index = %d\n", universal_index); fflush(stdout);
  
  // We can't test this besause 
  // assert(universal_index >=1 && universal_index <= MAX_RATES);

  // Can't assert this either because we might be doing this for another rate.
  // assert(available_rates[universal_index] == 0);

  available_rates[rate] = universal_rates[universal_index];

  if (dbg) printf("available_rates[%d] = %s\n", rate, available_rates[rate]);
}

int
ratestr_to_universal_index(char *ratestr) 
{
  int i = 1;
  // Find the equal sign (we skip ahead because we know it can't
  // be before 
  int equal_index = 11;
#if 0
  int end = equal_index + 3;
  for (; equal_index < end; i++) {
    if (config_str[i] == '=') {
      break;
    }
  }
  assert(equal_index != end);
#endif

  if (universal_rates[1] == NULL) {
    universal_rates_init();
  }

  // if (dbg) printf("equal_index = %d\n", equal_index);
  for (i = 1; i<=MAX_RATES; i++) {
    if (strncmp(ratestr, universal_rates[i], equal_index) == 0) {
      if (dbg) printf("found %s at %d\n", ratestr, i);
      return i;
    }
  }
  return -1;
}


void
universal_index_to_str(int i, char *ratestr)
{
  assert(i >= 1 && i <= MAX_RATES);

  // printf("universal_rates[1] = %p\n", universal_rates[1]); 
  if (universal_rates[1] == NULL) {
    universal_rates_init();
  }
  strcpy(ratestr, universal_rates[i]);
}

// TODO: This approach doesn't work if we want to use something
//       like 1 stream and 16 rates. We get 20/40 rates instead of SGI/LGI.

void
rateindex_to_str(int i, char *ratestr)
{
  strcpy(ratestr, available_rates[i]);
}

// TODO: This approach doesn't work if we want to use something
//       like 1 stream and 16 rates. We get 20/40 rates instead of SGI/LGI.

void
rateindex_to_str_new(int i, int streams, int widths, int rates, char *ratestr)
{
  int index = 0;
  int mcs = 0;
  int guard = 0;
  int width = 0;

  // int num_widths = sizeof(width_strs)/sizeof(char *);

  index = i-1;
  mcs = (index % (rates/(streams+1)));
  guard = (index / (rates/(streams+1))) % 2;
  width = (index / (rates/widths));


  if (dbg) printf("rates/streams = %d widths = %d index = %d mcs = %d guard = %d width = %d\n", (rates/streams), widths, index, mcs, guard, width);

  sprintf(ratestr, "%s-%s-%s=%s", mcs_strs[mcs], guard_strs[guard], width_strs[width], phyrate_strs[index]);

  fflush(stdout);
}

