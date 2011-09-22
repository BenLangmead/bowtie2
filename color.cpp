/*
 * Copyright 2011, Ben Langmead <blangmea@jhsph.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
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
