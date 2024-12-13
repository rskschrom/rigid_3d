#include <string>
#include "params.h"
#include "state.h"

void init_body(Params p, State &s, int nx, int ny);
void init_body_file(Params p, State &s, std::string fname);
