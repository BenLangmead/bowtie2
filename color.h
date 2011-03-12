/*
 * color.h
 *
 *  Created on: Oct 18, 2009
 *      Author: Ben Langmead
 */

#ifndef COLOR_H_
#define COLOR_H_

#include <string>

enum {
	COLOR_RED = 1,
	COLOR_GREEN,
	COLOR_YELLOW,
	COLOR_BLUE,
	COLOR_WHITE = 7
};


void appendConsoleColor(std::string& s, int color);
void setConsoleColor(int color);
void appendColor(std::string& s, char color);
void printColor(char color);

#endif /* COLOR_H_ */
