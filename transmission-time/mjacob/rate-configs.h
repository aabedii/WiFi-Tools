#include <assert.h>
#include <stdbool.h>
#include "rates.h"

char **read_rate_configs(char *filename, int *count);

int file_config_name_to_index(char *filename);

bool filename_in_rate_config(char *filename, char **rate_configs);

/* Given a rate string figure out the filename for that rate configuration */
int ratestr_to_filename(char *ratestr, char *filename);

void print_rate_info(int index);
