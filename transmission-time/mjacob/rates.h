#include <stdbool.h>

#define RATE_STR_MAX     (20)
#define MAX_RATES        (96)

// We have universal rates 1 ... MAX_RATES (i.e., 1...96)
// These stay the same not matter what rates are available in
// the configuration being used.

// Available rates is the subset of rates that is available
// in the current configuration of the system being used
// E.g., 1..16 might be 1 stream, with only 20 MHz (40 MHz not available)
// There is a mapping from an available rate = 1 to a univeral rate = x
// and from an available rate 1 to universal rate str = s.
// All strings are printed using the universal string

extern int RATES;

void rateindex_to_str(int i, char *ratestr);
void rateindex_to_str_new(int i, int streams, int widths, int rates, char *ratestr);

void available_rate_init(int rate, int mcs, int guard, int width);
void universal_index_to_str(int i, char *ratestr);
int ratestr_to_universal_index(char *config_str);

char *available_rate_str(int rate);

char *mcs_to_str(int mcs);
char *guard_to_str(int guard);
char *width_to_str(int width);

int mcs_str_to_index(char *mcs);
int guard_str_to_index(char *guard);
int width_str_to_index(char *width);


bool ratestr_is_in_available_rates(char *ratestr);
int ratestr_to_available_index(char *ratestr);

