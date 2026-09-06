/* Standalone oracle for tests/test_sequence.py.
 *
 * The bodies of method_composition and method_longest_orf from the original
 * Ruby transrate's C extension (ext/transrate/transrate.c), with the Ruby
 * VALUE wrappers removed and nothing else changed. The extension itself was
 * deleted when the port completed; this copy is kept deliberately, as the
 * reference the Python implementation is diffed against. Do not "tidy" it --
 * its value is that it is not our code.
 *
 * Reads one sequence per line on stdin, prints
 *   <a> <c> <g> <t> <n> <longest_orf>
 * per line. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int base_counts[5];
static int dibase_counts[25];

static void composition(const char *seq, int len) {
  int i, idx;
  char base, prevbase;
  for (i = 0; i < 5; i++) base_counts[i] = 0;
  for (i = 0; i < 25; i++) dibase_counts[i] = 0;
  for (i = 0; i < len; i++) {
    base = seq[i];
    if (base > 90) base -= 32;
    switch (base) {
      case 'A': idx = 0; break;
      case 'C': idx = 1; break;
      case 'G': idx = 2; break;
      case 'T': idx = 3; break;
      default:  idx = 4; break;
    }
    base_counts[idx]++;
    if (i > 0) {
      prevbase = seq[i-1];
      if (prevbase > 90) prevbase -= 32;
      switch (prevbase) {
        case 'A': idx = idx;      break;
        case 'C': idx = idx + 5;  break;
        case 'G': idx = idx + 10; break;
        case 'T': idx = idx + 15; break;
        default:  idx = idx + 20; break;
      }
      dibase_counts[idx]++;
    }
  }
}

static int longest_orf(const char *str, int sl) {
  int i, longest = 0;
  int len[3];
  for (i = 0; i < 3; i++) len[i] = 0;
  for (i = 0; i < sl - 2; i++) {
    if (str[i] == 'A' && str[i+1] == 'T' && str[i+2] == 'G') {
      if (len[i%3] >= 0) { len[i%3]++; } else { len[i%3] = 1; }
    } else {
      if (str[i] == 'T' &&
          ((str[i+1] == 'A' && str[i+2] == 'G') ||
           (str[i+1] == 'A' && str[i+2] == 'A') ||
           (str[i+1] == 'G' && str[i+2] == 'A'))) {
        if (len[i%3] > longest) longest = len[i%3];
        len[i%3] = -1;
      } else {
        if (len[i%3] >= 0) len[i%3]++;
      }
    }
  }
  for (i = 0; i < 3; i++) if (len[i%3] > longest) longest = len[i%3];

  for (i = 0; i < 3; i++) len[i] = 0;
  for (i = sl - 1; i >= 2; i--) {
    if (str[i] == 'T' && str[i-1] == 'A' && str[i-2] == 'C') {
      if (len[i%3] >= 0) { len[i%3]++; } else { len[i%3] = 1; }
    } else {
      if (str[i] == 'A' &&
          ((str[i-1] == 'T' && str[i-2] == 'C') ||
           (str[i-1] == 'T' && str[i-2] == 'T') ||
           (str[i-1] == 'C' && str[i-2] == 'T'))) {
        if (len[i%3] > longest) longest = len[i%3];
        len[i%3] = -1;
      } else {
        if (len[i%3] >= 0) len[i%3]++;
      }
    }
  }
  for (i = 0; i < 3; i++) if (len[i%3] > longest) longest = len[i%3];
  return longest;
}

int main(void) {
  char *line = NULL;
  size_t cap = 0;
  ssize_t got;
  while ((got = getline(&line, &cap, stdin)) != -1) {
    while (got > 0 && (line[got-1] == '\n' || line[got-1] == '\r')) line[--got] = '\0';
    composition(line, (int)got);
    printf("%d %d %d %d %d %d\n",
           base_counts[0], base_counts[1], base_counts[2],
           base_counts[3], base_counts[4], longest_orf(line, (int)got));
  }
  free(line);
  return 0;
}
