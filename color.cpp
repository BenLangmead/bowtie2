/*
 * color.cpp
 *
 *  Created on: Oct 18, 2009
 *      Author: Ben Langmead
 */

#include <iostream>
#include <sstream>
#include <string>
#include "color.h"

using namespace std;

/**
 * Set the console color.
 */
void setConsoleColor(int color) {
	cout << (char)0x1B << "[" << 0 << ";" << color + 30 << ";" << 0 + 40 << "m";
}

/**
 * Set the console color.
 */
void appendConsoleColor(string& s, int color) {
	s.push_back((char)0x1B);
	s.append("[0;");
	ostringstream ss;
	ss << (color + 30);
	s.append(ss.str());
	s.append(";40m");
}

/**
 * Print color character in the appropriate color to the console.
 * Console must support color.
 */
void printColor(char color) {
	char ch = ' ';
	switch(color) {
		case 'A':
		case '0':
		case 0:
			setConsoleColor(COLOR_BLUE);
			ch = '0'; break;
		case 'C':
		case '1':
		case 1:
			setConsoleColor(COLOR_GREEN);
			ch = '1'; break;
		case 'G':
		case '2':
		case 2:
			setConsoleColor(COLOR_YELLOW);
			ch = '2'; break;
		case 'T':
		case '3':
		case 3:
			setConsoleColor(COLOR_RED);
			ch = '3'; break;
		case 'N':
		case '4':
		case '.':
		case 4:
			setConsoleColor(COLOR_WHITE);
			ch = '.'; break;
		default:
			setConsoleColor(COLOR_WHITE);
			break;
	}
	cout << ch;
	setConsoleColor(COLOR_WHITE);
}

/**
 * Print color character in the appropriate color to the console.
 * Console must support color.
 */
void appendColor(string& s, char color) {
	char ch = ' ';
	switch(color) {
		case 'A':
		case '0':
		case 0:
			appendConsoleColor(s, COLOR_BLUE);
			ch = '0'; break;
		case 'C':
		case '1':
		case 1:
			appendConsoleColor(s, COLOR_GREEN);
			ch = '1'; break;
		case 'G':
		case '2':
		case 2:
			appendConsoleColor(s, COLOR_YELLOW);
			ch = '2'; break;
		case 'T':
		case '3':
		case 3:
			appendConsoleColor(s, COLOR_RED);
			ch = '3'; break;
		case 'N':
		case '4':
		case '.':
		case 4:
			appendConsoleColor(s, COLOR_WHITE);
			ch = '.'; break;
		default:
			appendConsoleColor(s, COLOR_WHITE);
			break;
	}
	s.push_back(ch);
	appendConsoleColor(s, COLOR_WHITE);
}
