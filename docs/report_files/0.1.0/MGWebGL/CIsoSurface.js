// File Name: CIsoSurface.cpp
// Last Modified: 16/2/2015
// Author: Raghavendra Chandrashekara (based on source code provided
// by Paul Bourke and Cory Gene Bloyd)
// Email: rc99@doc.ic.ac.uk, rchandrashekara@hotmail.com
//
// Description: This is the implementation file for the CIsoSurface class.
//
// Extra crystallograhic and export code added by Paul Emsley and Kevin Cowtan.
//
// Converted to JavaScript by Stuart McNicholas.

// Should require returnVertices, returnNormals_new and returnIndices.

var m_edgeTable = [
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
];

var m_triTable = [
	[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 1, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[1, 8, 3, 9, 8, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 8, 3, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[9, 2, 10, 0, 2, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[2, 8, 3, 2, 10, 8, 10, 9, 8, -1, -1, -1, -1, -1, -1, -1],
	[3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 11, 2, 8, 11, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[1, 9, 0, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[1, 11, 2, 1, 9, 11, 9, 8, 11, -1, -1, -1, -1, -1, -1, -1],
	[3, 10, 1, 11, 10, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 10, 1, 0, 8, 10, 8, 11, 10, -1, -1, -1, -1, -1, -1, -1],
	[3, 9, 0, 3, 11, 9, 11, 10, 9, -1, -1, -1, -1, -1, -1, -1],
	[9, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[4, 3, 0, 7, 3, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 1, 9, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[4, 1, 9, 4, 7, 1, 7, 3, 1, -1, -1, -1, -1, -1, -1, -1],
	[1, 2, 10, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[3, 4, 7, 3, 0, 4, 1, 2, 10, -1, -1, -1, -1, -1, -1, -1],
	[9, 2, 10, 9, 0, 2, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
	[2, 10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4, -1, -1, -1, -1],
	[8, 4, 7, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[11, 4, 7, 11, 2, 4, 2, 0, 4, -1, -1, -1, -1, -1, -1, -1],
	[9, 0, 1, 8, 4, 7, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
	[4, 7, 11, 9, 4, 11, 9, 11, 2, 9, 2, 1, -1, -1, -1, -1],
	[3, 10, 1, 3, 11, 10, 7, 8, 4, -1, -1, -1, -1, -1, -1, -1],
	[1, 11, 10, 1, 4, 11, 1, 0, 4, 7, 11, 4, -1, -1, -1, -1],
	[4, 7, 8, 9, 0, 11, 9, 11, 10, 11, 0, 3, -1, -1, -1, -1],
	[4, 7, 11, 4, 11, 9, 9, 11, 10, -1, -1, -1, -1, -1, -1, -1],
	[9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[9, 5, 4, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 5, 4, 1, 5, 0, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[8, 5, 4, 8, 3, 5, 3, 1, 5, -1, -1, -1, -1, -1, -1, -1],
	[1, 2, 10, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[3, 0, 8, 1, 2, 10, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
	[5, 2, 10, 5, 4, 2, 4, 0, 2, -1, -1, -1, -1, -1, -1, -1],
	[2, 10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8, -1, -1, -1, -1],
	[9, 5, 4, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 11, 2, 0, 8, 11, 4, 9, 5, -1, -1, -1, -1, -1, -1, -1],
	[0, 5, 4, 0, 1, 5, 2, 3, 11, -1, -1, -1, -1, -1, -1, -1],
	[2, 1, 5, 2, 5, 8, 2, 8, 11, 4, 8, 5, -1, -1, -1, -1],
	[10, 3, 11, 10, 1, 3, 9, 5, 4, -1, -1, -1, -1, -1, -1, -1],
	[4, 9, 5, 0, 8, 1, 8, 10, 1, 8, 11, 10, -1, -1, -1, -1],
	[5, 4, 0, 5, 0, 11, 5, 11, 10, 11, 0, 3, -1, -1, -1, -1],
	[5, 4, 8, 5, 8, 10, 10, 8, 11, -1, -1, -1, -1, -1, -1, -1],
	[9, 7, 8, 5, 7, 9, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[9, 3, 0, 9, 5, 3, 5, 7, 3, -1, -1, -1, -1, -1, -1, -1],
	[0, 7, 8, 0, 1, 7, 1, 5, 7, -1, -1, -1, -1, -1, -1, -1],
	[1, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[9, 7, 8, 9, 5, 7, 10, 1, 2, -1, -1, -1, -1, -1, -1, -1],
	[10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3, -1, -1, -1, -1],
	[8, 0, 2, 8, 2, 5, 8, 5, 7, 10, 5, 2, -1, -1, -1, -1],
	[2, 10, 5, 2, 5, 3, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1],
	[7, 9, 5, 7, 8, 9, 3, 11, 2, -1, -1, -1, -1, -1, -1, -1],
	[9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7, 11, -1, -1, -1, -1],
	[2, 3, 11, 0, 1, 8, 1, 7, 8, 1, 5, 7, -1, -1, -1, -1],
	[11, 2, 1, 11, 1, 7, 7, 1, 5, -1, -1, -1, -1, -1, -1, -1],
	[9, 5, 8, 8, 5, 7, 10, 1, 3, 10, 3, 11, -1, -1, -1, -1],
	[5, 7, 0, 5, 0, 9, 7, 11, 0, 1, 0, 10, 11, 10, 0, -1],
	[11, 10, 0, 11, 0, 3, 10, 5, 0, 8, 0, 7, 5, 7, 0, -1],
	[11, 10, 5, 7, 11, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 8, 3, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[9, 0, 1, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[1, 8, 3, 1, 9, 8, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
	[1, 6, 5, 2, 6, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[1, 6, 5, 1, 2, 6, 3, 0, 8, -1, -1, -1, -1, -1, -1, -1],
	[9, 6, 5, 9, 0, 6, 0, 2, 6, -1, -1, -1, -1, -1, -1, -1],
	[5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8, -1, -1, -1, -1],
	[2, 3, 11, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[11, 0, 8, 11, 2, 0, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
	[0, 1, 9, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1, -1, -1, -1],
	[5, 10, 6, 1, 9, 2, 9, 11, 2, 9, 8, 11, -1, -1, -1, -1],
	[6, 3, 11, 6, 5, 3, 5, 1, 3, -1, -1, -1, -1, -1, -1, -1],
	[0, 8, 11, 0, 11, 5, 0, 5, 1, 5, 11, 6, -1, -1, -1, -1],
	[3, 11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9, -1, -1, -1, -1],
	[6, 5, 9, 6, 9, 11, 11, 9, 8, -1, -1, -1, -1, -1, -1, -1],
	[5, 10, 6, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[4, 3, 0, 4, 7, 3, 6, 5, 10, -1, -1, -1, -1, -1, -1, -1],
	[1, 9, 0, 5, 10, 6, 8, 4, 7, -1, -1, -1, -1, -1, -1, -1],
	[10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4, -1, -1, -1, -1],
	[6, 1, 2, 6, 5, 1, 4, 7, 8, -1, -1, -1, -1, -1, -1, -1],
	[1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7, -1, -1, -1, -1],
	[8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6, -1, -1, -1, -1],
	[7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9, -1],
	[3, 11, 2, 7, 8, 4, 10, 6, 5, -1, -1, -1, -1, -1, -1, -1],
	[5, 10, 6, 4, 7, 2, 4, 2, 0, 2, 7, 11, -1, -1, -1, -1],
	[0, 1, 9, 4, 7, 8, 2, 3, 11, 5, 10, 6, -1, -1, -1, -1],
	[9, 2, 1, 9, 11, 2, 9, 4, 11, 7, 11, 4, 5, 10, 6, -1],
	[8, 4, 7, 3, 11, 5, 3, 5, 1, 5, 11, 6, -1, -1, -1, -1],
	[5, 1, 11, 5, 11, 6, 1, 0, 11, 7, 11, 4, 0, 4, 11, -1],
	[0, 5, 9, 0, 6, 5, 0, 3, 6, 11, 6, 3, 8, 4, 7, -1],
	[6, 5, 9, 6, 9, 11, 4, 7, 9, 7, 11, 9, -1, -1, -1, -1],
	[10, 4, 9, 6, 4, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[4, 10, 6, 4, 9, 10, 0, 8, 3, -1, -1, -1, -1, -1, -1, -1],
	[10, 0, 1, 10, 6, 0, 6, 4, 0, -1, -1, -1, -1, -1, -1, -1],
	[8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1, 10, -1, -1, -1, -1],
	[1, 4, 9, 1, 2, 4, 2, 6, 4, -1, -1, -1, -1, -1, -1, -1],
	[3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4, -1, -1, -1, -1],
	[0, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[8, 3, 2, 8, 2, 4, 4, 2, 6, -1, -1, -1, -1, -1, -1, -1],
	[10, 4, 9, 10, 6, 4, 11, 2, 3, -1, -1, -1, -1, -1, -1, -1],
	[0, 8, 2, 2, 8, 11, 4, 9, 10, 4, 10, 6, -1, -1, -1, -1],
	[3, 11, 2, 0, 1, 6, 0, 6, 4, 6, 1, 10, -1, -1, -1, -1],
	[6, 4, 1, 6, 1, 10, 4, 8, 1, 2, 1, 11, 8, 11, 1, -1],
	[9, 6, 4, 9, 3, 6, 9, 1, 3, 11, 6, 3, -1, -1, -1, -1],
	[8, 11, 1, 8, 1, 0, 11, 6, 1, 9, 1, 4, 6, 4, 1, -1],
	[3, 11, 6, 3, 6, 0, 0, 6, 4, -1, -1, -1, -1, -1, -1, -1],
	[6, 4, 8, 11, 6, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[7, 10, 6, 7, 8, 10, 8, 9, 10, -1, -1, -1, -1, -1, -1, -1],
	[0, 7, 3, 0, 10, 7, 0, 9, 10, 6, 7, 10, -1, -1, -1, -1],
	[10, 6, 7, 1, 10, 7, 1, 7, 8, 1, 8, 0, -1, -1, -1, -1],
	[10, 6, 7, 10, 7, 1, 1, 7, 3, -1, -1, -1, -1, -1, -1, -1],
	[1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7, -1, -1, -1, -1],
	[2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9, -1],
	[7, 8, 0, 7, 0, 6, 6, 0, 2, -1, -1, -1, -1, -1, -1, -1],
	[7, 3, 2, 6, 7, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[2, 3, 11, 10, 6, 8, 10, 8, 9, 8, 6, 7, -1, -1, -1, -1],
	[2, 0, 7, 2, 7, 11, 0, 9, 7, 6, 7, 10, 9, 10, 7, -1],
	[1, 8, 0, 1, 7, 8, 1, 10, 7, 6, 7, 10, 2, 3, 11, -1],
	[11, 2, 1, 11, 1, 7, 10, 6, 1, 6, 7, 1, -1, -1, -1, -1],
	[8, 9, 6, 8, 6, 7, 9, 1, 6, 11, 6, 3, 1, 3, 6, -1],
	[0, 9, 1, 11, 6, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[7, 8, 0, 7, 0, 6, 3, 11, 0, 11, 6, 0, -1, -1, -1, -1],
	[7, 11, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[3, 0, 8, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 1, 9, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[8, 1, 9, 8, 3, 1, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
	[10, 1, 2, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[1, 2, 10, 3, 0, 8, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
	[2, 9, 0, 2, 10, 9, 6, 11, 7, -1, -1, -1, -1, -1, -1, -1],
	[6, 11, 7, 2, 10, 3, 10, 8, 3, 10, 9, 8, -1, -1, -1, -1],
	[7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1, -1, -1, -1, -1],
	[2, 7, 6, 2, 3, 7, 0, 1, 9, -1, -1, -1, -1, -1, -1, -1],
	[1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6, -1, -1, -1, -1],
	[10, 7, 6, 10, 1, 7, 1, 3, 7, -1, -1, -1, -1, -1, -1, -1],
	[10, 7, 6, 1, 7, 10, 1, 8, 7, 1, 0, 8, -1, -1, -1, -1],
	[0, 3, 7, 0, 7, 10, 0, 10, 9, 6, 10, 7, -1, -1, -1, -1],
	[7, 6, 10, 7, 10, 8, 8, 10, 9, -1, -1, -1, -1, -1, -1, -1],
	[6, 8, 4, 11, 8, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[3, 6, 11, 3, 0, 6, 0, 4, 6, -1, -1, -1, -1, -1, -1, -1],
	[8, 6, 11, 8, 4, 6, 9, 0, 1, -1, -1, -1, -1, -1, -1, -1],
	[9, 4, 6, 9, 6, 3, 9, 3, 1, 11, 3, 6, -1, -1, -1, -1],
	[6, 8, 4, 6, 11, 8, 2, 10, 1, -1, -1, -1, -1, -1, -1, -1],
	[1, 2, 10, 3, 0, 11, 0, 6, 11, 0, 4, 6, -1, -1, -1, -1],
	[4, 11, 8, 4, 6, 11, 0, 2, 9, 2, 10, 9, -1, -1, -1, -1],
	[10, 9, 3, 10, 3, 2, 9, 4, 3, 11, 3, 6, 4, 6, 3, -1],
	[8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1],
	[0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8, -1, -1, -1, -1],
	[1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1],
	[8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 10, 1, -1, -1, -1, -1],
	[10, 1, 0, 10, 0, 6, 6, 0, 4, -1, -1, -1, -1, -1, -1, -1],
	[4, 6, 3, 4, 3, 8, 6, 10, 3, 0, 3, 9, 10, 9, 3, -1],
	[10, 9, 4, 6, 10, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[4, 9, 5, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 8, 3, 4, 9, 5, 11, 7, 6, -1, -1, -1, -1, -1, -1, -1],
	[5, 0, 1, 5, 4, 0, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
	[11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5, -1, -1, -1, -1],
	[9, 5, 4, 10, 1, 2, 7, 6, 11, -1, -1, -1, -1, -1, -1, -1],
	[6, 11, 7, 1, 2, 10, 0, 8, 3, 4, 9, 5, -1, -1, -1, -1],
	[7, 6, 11, 5, 4, 10, 4, 2, 10, 4, 0, 2, -1, -1, -1, -1],
	[3, 4, 8, 3, 5, 4, 3, 2, 5, 10, 5, 2, 11, 7, 6, -1],
	[7, 2, 3, 7, 6, 2, 5, 4, 9, -1, -1, -1, -1, -1, -1, -1],
	[9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7, -1, -1, -1, -1],
	[3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0, -1, -1, -1, -1],
	[6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8, -1],
	[9, 5, 4, 10, 1, 6, 1, 7, 6, 1, 3, 7, -1, -1, -1, -1],
	[1, 6, 10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4, -1],
	[4, 0, 10, 4, 10, 5, 0, 3, 10, 6, 10, 7, 3, 7, 10, -1],
	[7, 6, 10, 7, 10, 8, 5, 4, 10, 4, 8, 10, -1, -1, -1, -1],
	[6, 9, 5, 6, 11, 9, 11, 8, 9, -1, -1, -1, -1, -1, -1, -1],
	[3, 6, 11, 0, 6, 3, 0, 5, 6, 0, 9, 5, -1, -1, -1, -1],
	[0, 11, 8, 0, 5, 11, 0, 1, 5, 5, 6, 11, -1, -1, -1, -1],
	[6, 11, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1, -1, -1, -1, -1],
	[1, 2, 10, 9, 5, 11, 9, 11, 8, 11, 5, 6, -1, -1, -1, -1],
	[0, 11, 3, 0, 6, 11, 0, 9, 6, 5, 6, 9, 1, 2, 10, -1],
	[11, 8, 5, 11, 5, 6, 8, 0, 5, 10, 5, 2, 0, 2, 5, -1],
	[6, 11, 3, 6, 3, 5, 2, 10, 3, 10, 5, 3, -1, -1, -1, -1],
	[5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2, -1, -1, -1, -1],
	[9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1, -1, -1, -1, -1],
	[1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8, -1],
	[1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[1, 3, 6, 1, 6, 10, 3, 8, 6, 5, 6, 9, 8, 9, 6, -1],
	[10, 1, 0, 10, 0, 6, 9, 5, 0, 5, 6, 0, -1, -1, -1, -1],
	[0, 3, 8, 5, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[10, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[11, 5, 10, 7, 5, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[11, 5, 10, 11, 7, 5, 8, 3, 0, -1, -1, -1, -1, -1, -1, -1],
	[5, 11, 7, 5, 10, 11, 1, 9, 0, -1, -1, -1, -1, -1, -1, -1],
	[10, 7, 5, 10, 11, 7, 9, 8, 1, 8, 3, 1, -1, -1, -1, -1],
	[11, 1, 2, 11, 7, 1, 7, 5, 1, -1, -1, -1, -1, -1, -1, -1],
	[0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2, 11, -1, -1, -1, -1],
	[9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 11, 7, -1, -1, -1, -1],
	[7, 5, 2, 7, 2, 11, 5, 9, 2, 3, 2, 8, 9, 8, 2, -1],
	[2, 5, 10, 2, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1],
	[8, 2, 0, 8, 5, 2, 8, 7, 5, 10, 2, 5, -1, -1, -1, -1],
	[9, 0, 1, 5, 10, 3, 5, 3, 7, 3, 10, 2, -1, -1, -1, -1],
	[9, 8, 2, 9, 2, 1, 8, 7, 2, 10, 2, 5, 7, 5, 2, -1],
	[1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1, -1, -1, -1, -1],
	[9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1, -1, -1, -1, -1],
	[9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[5, 8, 4, 5, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1],
	[5, 0, 4, 5, 11, 0, 5, 10, 11, 11, 3, 0, -1, -1, -1, -1],
	[0, 1, 9, 8, 4, 10, 8, 10, 11, 10, 4, 5, -1, -1, -1, -1],
	[10, 11, 4, 10, 4, 5, 11, 3, 4, 9, 4, 1, 3, 1, 4, -1],
	[2, 5, 1, 2, 8, 5, 2, 11, 8, 4, 5, 8, -1, -1, -1, -1],
	[0, 4, 11, 0, 11, 3, 4, 5, 11, 2, 11, 1, 5, 1, 11, -1],
	[0, 2, 5, 0, 5, 9, 2, 11, 5, 4, 5, 8, 11, 8, 5, -1],
	[9, 4, 5, 2, 11, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[2, 5, 10, 3, 5, 2, 3, 4, 5, 3, 8, 4, -1, -1, -1, -1],
	[5, 10, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1, -1, -1, -1, -1],
	[3, 10, 2, 3, 5, 10, 3, 8, 5, 4, 5, 8, 0, 1, 9, -1],
	[5, 10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2, -1, -1, -1, -1],
	[8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1, -1, -1, -1, -1],
	[0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5, -1, -1, -1, -1],
	[9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[4, 11, 7, 4, 9, 11, 9, 10, 11, -1, -1, -1, -1, -1, -1, -1],
	[0, 8, 3, 4, 9, 7, 9, 11, 7, 9, 10, 11, -1, -1, -1, -1],
	[1, 10, 11, 1, 11, 4, 1, 4, 0, 7, 4, 11, -1, -1, -1, -1],
	[3, 1, 4, 3, 4, 8, 1, 10, 4, 7, 4, 11, 10, 11, 4, -1],
	[4, 11, 7, 9, 11, 4, 9, 2, 11, 9, 1, 2, -1, -1, -1, -1],
	[9, 7, 4, 9, 11, 7, 9, 1, 11, 2, 11, 1, 0, 8, 3, -1],
	[11, 7, 4, 11, 4, 2, 2, 4, 0, -1, -1, -1, -1, -1, -1, -1],
	[11, 7, 4, 11, 4, 2, 8, 3, 4, 3, 2, 4, -1, -1, -1, -1],
	[2, 9, 10, 2, 7, 9, 2, 3, 7, 7, 4, 9, -1, -1, -1, -1],
	[9, 10, 7, 9, 7, 4, 10, 2, 7, 8, 7, 0, 2, 0, 7, -1],
	[3, 7, 10, 3, 10, 2, 7, 4, 10, 1, 10, 0, 4, 0, 10, -1],
	[1, 10, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1, -1, -1, -1, -1],
	[4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1, -1, -1, -1, -1],
	[4, 0, 3, 7, 4, 3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[4, 8, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[9, 10, 8, 10, 11, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[3, 0, 9, 3, 9, 11, 11, 9, 10, -1, -1, -1, -1, -1, -1, -1],
	[0, 1, 10, 0, 10, 8, 8, 10, 11, -1, -1, -1, -1, -1, -1, -1],
	[3, 1, 10, 11, 3, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[1, 2, 11, 1, 11, 9, 9, 11, 8, -1, -1, -1, -1, -1, -1, -1],
	[3, 0, 9, 3, 9, 11, 1, 2, 9, 2, 11, 9, -1, -1, -1, -1],
	[0, 2, 11, 8, 0, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[3, 2, 11, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[2, 3, 8, 2, 8, 10, 10, 8, 9, -1, -1, -1, -1, -1, -1, -1],
	[9, 10, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[2, 3, 8, 2, 8, 10, 0, 1, 8, 1, 10, 8, -1, -1, -1, -1],
	[1, 10, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[1, 3, 8, 9, 1, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
	[-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1]
];

function CIsoSurface() {
    this.m_fCellLengthX = 0;
    this.m_fCellLengthY = 0;
    this.m_fCellLengthZ = 0;
    this.m_nCellsX = 0;
    this.m_nCellsY = 0;
    this.m_nCellsZ = 0;
    this.m_nTriangles = 0;
    this.m_nTriangles_clipped = 0;
    this.m_nNormals = 0;
    this.m_nVertices = 0;
    this.m_ppt3dVertices = [];
    this.m_piTriangleIndices = [];
    this.m_pvec3dNormals = [];
    this.m_ptScalarField = [];
    this.m_tIsoLevel = 0;
    this.m_bValidSurface = false;
    this.m_i2pt3idVertices = {};
    this.m_trivecTriangles = [];
    this.chicken_vertices = [];
    this.cell = null;
}

function VectorDiff(p0,p1){
    var pnew = [p0[0]-p1[0],p0[1]-p1[1],p0[2]-p1[2],1.0];
    return pnew;
}

function VectorSum(p0,p1){
    var pnew = [p0[0]+p1[0],p0[1]+p1[1],p0[2]+p1[2],1.0];
    return pnew;
}

function VectorNormalized(p){
    var fac = Math.sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    var pnew = [p[0]/fac,p[1]/fac,p[2]/fac];
    return pnew;
}

function PlaneNormalized(p){
    var fac = Math.sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    var pnew = [p[0]/fac,p[1]/fac,p[2]/fac,p[3]/fac];
    return pnew;
}

function VectorLength(p){
    return  Math.sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}

function printVec4(p){
    console.log(p[0]+" "+p[1]+" "+p[2]+" "+p[3]);
}

function printVec(p){
    console.log(p[0]+" "+p[1]+" "+p[2]);
}

function Quat(x,y,z,wi,angle_in){
    
    var q;
    if(wi==0){
        var xQ = Quat(1.0,0.0,0.0,1,x);
        var yQ = Quat(0.0,0.1,0.0,1,y);
        var zQ = Quat(0.0,0.0,1.0,1,z);
        q = quat4.create(xQ);
        quat4.multiply(xQ,yQ);
        quat4.multiply(xQ,zQ);
    }else{
        var angle = angle_in * Math.PI/180.0;
        var radius = Math.sqrt(x*x+ y*y+z*z);
        if(radius < 0.000000001){
            console.log("zero length in Quat\n");
            return;
        }
        var x2 = x / radius;
        var y2 = y / radius;
        var z2 = z / radius;

        var dval = [];
        dval[3] = Math.cos(angle/2.0);
        dval[0] = x2*Math.sin(angle/2.0);
        dval[1] = y2*Math.sin(angle/2.0);
        dval[2] = z2*Math.sin(angle/2.0);
        q = quat4.create(dval);
        
    }
    return q;
}

function Angle(A, B, C){

  var BA = VectorDiff(B,A);
  var CA = VectorDiff(C,A);
  var CB = VectorDiff(C,B);

  var ab = VectorLength(BA);
  var ac = VectorLength(CA);
  var bc = VectorLength(CB);

  var absq = ab*ab;
  var acsq = ac*ac;
  var bcsq = bc*bc;

  return  Math.acos((bcsq + absq - acsq)/(2*bc*ab));
}

function FloatVectorMult(f,p){
    var vec = [0.0,0.0,0.0,1.0];
    vec[0] += f*p[0];
    vec[1] += f*p[1];
    vec[2] += f*p[2];
    return vec;
}

function VectorNegate(p){
    var vec = [0.0,0.0,0.0,1.0];
    vec[0] -= p[0];
    vec[1] -= p[1];
    vec[2] -= p[2];
    return vec;
}

function VectorCopy(p){
    var vec = [0.0,0.0,0.0,1.0];
    vec[0] += p[0];
    vec[1] += p[1];
    vec[2] += p[2];
    return vec;
}

function PlaneGetNormal(p){
    var vec = [0.0,0.0,0.0,1.0];
    vec[0] += p[0];
    vec[1] += p[1];
    vec[2] += p[2];
    return vec;
}

function RotMxV(objrotmat, prim){
  var result = [0.0,0.0,0.0,1.0];
  var input = [prim[0], prim[1], prim[2], 1.0];

  for(var i=0;i<4;i++){
    result[i] = 0.0;
    for(var j=0;j<4;j++){
      result[i] += input[j]*objrotmat[i*4+j];
    }
  }
  
  return result;
}

function line_plane_intersection(plane,p0,p1){

    var intersec = [];

    var p = [p1[0]-p0[0],p1[1]-p0[1],p1[2]-p0[2]];

    var num = plane[0]*p0[0] + plane[1]*p0[1] + plane[2]*p0[2] + plane[3];

    var den = plane[0]*p[0] + plane[1]*p[1] + plane[2]*p[2];

    if(Math.abs(den) < 0.000001){
        return [1e+7,1e+7,1e+7,1e+7];// This is dodgy, better way?
    }
    var u = -num/den;

    intersec[0] = p0[0] + u*(p1[0]-p0[0]);
    intersec[1] = p0[1] + u*(p1[1]-p0[1]);
    intersec[2] = p0[2] + u*(p1[2]-p0[2]);
    intersec[3] = 1.0;

    return intersec;
}

function DotProduct(v1, v2){
  var v;

  v  = v1[0] * v2[0];
  v += v1[1] * v2[1];
  v += v1[2] * v2[2];
  return v;
}

function CrossProduct(v1, v2){
  var v = [0.0,0.0,0.0,1.0];
  v[0] = v1[1] * v2[2] - v2[1] * v1[2];
  v[1] = v1[2] * v2[0] - v2[2] * v1[0];
  v[2] = v1[0] * v2[1] - v2[0] * v1[1];
  return v;
}

function VectorMidPoint(v1, v2){
  var v = [0.0,0.0,0.0,1.0];
  v[0] = v1[0]+(v2[0]-v1[0])/2;
  v[1] = v1[1]+(v2[1]-v1[1])/2;
  v[2] = v1[2]+(v2[2]-v1[2])/2;
  return v;
}

function VectorArrayMidPoint(v1){
  var v = [0.0,0.0,0.0,1.0];
  if(v1.length==0) return v;
  for(var i=0;i<v1.length;i++){
    v = VectorSum(v,v1[i]);
  }
  v[0] /= v1.length;
  v[1] /= v1.length;
  v[2] /= v1.length;
  return v;
}

CIsoSurface.prototype.ScanAdjacenciesYZ = function(x, y, z, way, depth){

    var vert_stack  = [];
    var horiz_stack = [];

    var depth_orig = JSON.parse(JSON.stringify(depth));
    if(depth>200){
        //console.log("Warning, returning because of potentially bad deep recursion");
        return;
    }
    depth++;

    var px = x*this.m_fCellLengthX;
    var p1y, p1z, p2y, p2z;
    var val1, val2;

    var vert = 0;
    var horiz = 0;
    var drawn = 0;

    var nPointsInXDirection = (this.m_nCellsX + 1);
    var nPointsInSlice = nPointsInXDirection*(this.m_nCellsY + 1);

    var intersection;

    if(way==1){
        horiz = 1;
    } else {
        vert = 1;
    }

    if(y>0&&z>0&&y<this.m_nCellsZ-1&&z<this.m_nCellsZ-1){

        if(vert){
            if(this.vert_edge[z][y+1]==1)
                vert_stack.push([z,y+1]);
            if(this.vert_edge[z][y-1]==1)
                vert_stack.push([z,y-1]);
            if(this.horiz_edge[z][y]==1)
                horiz_stack.push([z,y]);
            if(this.horiz_edge[z+1][y]==1)
                horiz_stack.push([z+1,y]);
            if(this.horiz_edge[z][y-1]==1)
                horiz_stack.push([z,y-1]);
            if(this.horiz_edge[z+1][y-1]==1)
                horiz_stack.push([z+1,y-1]);
        }else{
            if(this.horiz_edge[z+1][y]==1)
                horiz_stack.push([z+1,y]);
            if(this.horiz_edge[z-1][y]==1)
                horiz_stack.push([z-1,y]);
            if(this.vert_edge[z][y]==1)
                vert_stack.push([z,y]);
            if(this.vert_edge[z-1][y]==1)
                vert_stack.push([z-1,y]);
            if(this.vert_edge[z][y+1]==1)
                vert_stack.push([z,y+1]);
            if(this.vert_edge[z-1][y+1]==1)
                vert_stack.push([z-1,y+1]);
        }


        if(vert_stack.length>0){
            var vert_stack_idx = vert_stack.length - 1;

            p1y = y*this.m_fCellLengthY;
            p1z = z*this.m_fCellLengthZ;
            val1 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x];
            if(vert){
                p2y = y*this.m_fCellLengthY;
                p2z = (z+1)*this.m_fCellLengthZ;
                val2 = this.m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + x];
            }else{
                p2y = (y+1)*this.m_fCellLengthY;
                p2z = z*this.m_fCellLengthZ;
                val2 = this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x];
            }

            var mu;
            mu = (this.m_tIsoLevel - val1)/(val2 - val1);
            var i1 = p1y + mu*(p2y - p1y);
            var i2 = p1z + mu*(p2z - p1z);

            while(vert_stack_idx>=0){
                this.chicken_vertices.push(px);
                this.chicken_vertices.push(i1);
                this.chicken_vertices.push(i2);

                p1y = vert_stack[vert_stack_idx][1]*this.m_fCellLengthY;
                p1z = vert_stack[vert_stack_idx][0]*this.m_fCellLengthZ;
                p2y = vert_stack[vert_stack_idx][1]*this.m_fCellLengthY;
                p2z = (vert_stack[vert_stack_idx][0]+1)*this.m_fCellLengthZ;

                val1 = this.m_ptScalarField[vert_stack[vert_stack_idx][0]*nPointsInSlice + vert_stack[vert_stack_idx][1]*nPointsInXDirection + x];
                val2 = this.m_ptScalarField[(vert_stack[vert_stack_idx][0]+1)*nPointsInSlice + vert_stack[vert_stack_idx][1]*nPointsInXDirection + x];

                var mu2;
                mu2 = (this.m_tIsoLevel - val1)/(val2 - val1);
                var i1_2 = p1y + mu2*(p2y - p1y);
                var i2_2 = p1z + mu2*(p2z - p1z);

                this.chicken_vertices.push(px);
                this.chicken_vertices.push(i1_2);
                this.chicken_vertices.push(i2_2);

                if(this.vert_edge[vert_stack[vert_stack_idx][0]][vert_stack[vert_stack_idx][1]]==1){
                    this.vert_edge[vert_stack[vert_stack_idx][0]][vert_stack[vert_stack_idx][1]]=0;
                    this.ScanAdjacenciesYZ(x,vert_stack[vert_stack_idx][1],vert_stack[vert_stack_idx][0],0,depth);
                    depth = JSON.parse(JSON.stringify(depth_orig));
                    //depth = depth_orig;
                }
                //vert_stack.pop_back();
                vert_stack_idx--;
            }
        }
        if(horiz_stack.length>0){

            var horiz_stack_idx = horiz_stack.length - 1;
            p1y = y*this.m_fCellLengthY;
            p1z = z*this.m_fCellLengthZ;
            val1 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x];
            if(vert){
                p2y = y*this.m_fCellLengthY;
                p2z = (z+1)*this.m_fCellLengthZ;
                val2 = this.m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + x];
            }else{
                p2y = (y+1)*this.m_fCellLengthY;
                p2z = z*this.m_fCellLengthZ;
                val2 = this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x];
            }

            var mu;
            mu = (this.m_tIsoLevel - val1)/(val2 - val1);
            var i1 = p1y + mu*(p2y - p1y);
            var i2 = p1z + mu*(p2z - p1z);

            while(horiz_stack_idx>=0){
                this.chicken_vertices.push(px);
                this.chicken_vertices.push(i1);
                this.chicken_vertices.push(i2);

                p1y = horiz_stack[horiz_stack_idx][1]*this.m_fCellLengthY;
                p1z = horiz_stack[horiz_stack_idx][0]*this.m_fCellLengthZ;
                p2y = (horiz_stack[horiz_stack_idx][1]+1)*this.m_fCellLengthY;
                p2z = horiz_stack[horiz_stack_idx][0]*this.m_fCellLengthZ;

                val1 = this.m_ptScalarField[horiz_stack[horiz_stack_idx][0]*nPointsInSlice + horiz_stack[horiz_stack_idx][1]*nPointsInXDirection + x];
                val2 = this.m_ptScalarField[horiz_stack[horiz_stack_idx][0]*nPointsInSlice + (horiz_stack[horiz_stack_idx][1]+1)*nPointsInXDirection + x];

                var mu2;
                mu2 = (this.m_tIsoLevel - val1)/(val2 - val1);
                var i1_2 = p1y + mu2*(p2y - p1y);
                var i2_2 = p1z + mu2*(p2z - p1z);

                this.chicken_vertices.push(px);
                this.chicken_vertices.push(i1_2);
                this.chicken_vertices.push(i2_2);

                if(this.horiz_edge[horiz_stack[horiz_stack_idx][0]][horiz_stack[horiz_stack_idx][1]]==1){
                    this.horiz_edge[horiz_stack[horiz_stack_idx][0]][horiz_stack[horiz_stack_idx][1]]=0;
                    this.ScanAdjacenciesYZ(x,horiz_stack[horiz_stack_idx][1],horiz_stack[horiz_stack_idx][0],1,depth);
                    depth = JSON.parse(JSON.stringify(depth_orig));
                    //depth = depth_orig;
                }
                //horiz_stack.pop_back();
                horiz_stack_idx--;
            }
        }
    }
}

CIsoSurface.prototype.ScanAdjacenciesXZ = function(x, y, z, way, depth){

    var vert_stack  = [];
    var horiz_stack = [];

    var depth_orig = depth;
    if(depth>200){
        //console.log("Warning, returning because of potentially bad deep recursion");
        return;
    }
    depth++;

    var py = y*this.m_fCellLengthY;
    var p1x, p1z, p2x, p2z;
    var val1, val2;

    var vert = 0;
    var horiz = 0;
    var drawn = 0;

    var nPointsInXDirection = (this.m_nCellsX + 1);
    var nPointsInSlice = nPointsInXDirection*(this.m_nCellsY + 1);

    var intersection;

    if(way==1){
        horiz = 1;
    } else {
        vert = 1;
    }

    if(x>0&&z>0&&x<this.m_nCellsX-1&&z<this.m_nCellsZ-1){

        if(vert){
            if(this.vert_edge[z][x+1]==1)
                vert_stack.push([z,x+1]);
            if(this.vert_edge[z][x-1]==1)
                vert_stack.push([z,x-1]);
            if(this.horiz_edge[z][x]==1)
                horiz_stack.push([z,x]);
            if(this.horiz_edge[z+1][x]==1)
                horiz_stack.push([z+1,x]);
            if(this.horiz_edge[z][x-1]==1)
                horiz_stack.push([z,x-1]);
            if(this.horiz_edge[z+1][x-1]==1)
                horiz_stack.push([z+1,x-1]);
        }else{
            if(this.horiz_edge[z+1][x]==1)
                horiz_stack.push([z+1,x]);
            if(this.horiz_edge[z-1][x]==1)
                horiz_stack.push([z-1,x]);
            if(this.vert_edge[z][x]==1)
                vert_stack.push([z,x]);
            if(this.vert_edge[z-1][x]==1)
                vert_stack.push([z-1,x]);
            if(this.vert_edge[z][x+1]==1)
                vert_stack.push([z,x+1]);
            if(this.vert_edge[z-1][x+1]==1)
                vert_stack.push([z-1,x+1]);
        }


        if(vert_stack.length>0){

            var vert_stack_idx = vert_stack.length - 1;
            p1x = x*this.m_fCellLengthX;
            p1z = z*this.m_fCellLengthZ;
            val1 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x];
            if(vert){
                p2x = x*this.m_fCellLengthX;
                p2z = (z+1)*this.m_fCellLengthZ;
                val2 = this.m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + x];
            }else{
                p2x = (x+1)*this.m_fCellLengthX;
                p2z = z*this.m_fCellLengthZ;
                val2 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x+1];
            }

            var mu;
            mu = (this.m_tIsoLevel - val1)/(val2 - val1);
            var i0 = p1x + mu*(p2x - p1x);
            var i2 = p1z + mu*(p2z - p1z);

            while(vert_stack_idx>=0){
                this.chicken_vertices.push(i0);
                this.chicken_vertices.push(py);
                this.chicken_vertices.push(i2);

                p1x = vert_stack[vert_stack_idx][1]*this.m_fCellLengthX;
                p1z = vert_stack[vert_stack_idx][0]*this.m_fCellLengthZ;
                p2x = vert_stack[vert_stack_idx][1]*this.m_fCellLengthX;
                p2z = (vert_stack[vert_stack_idx][0]+1)*this.m_fCellLengthZ;

                val1 = this.m_ptScalarField[vert_stack[vert_stack_idx][0]*nPointsInSlice + y*nPointsInXDirection + vert_stack[vert_stack_idx][1]];
                val2 = this.m_ptScalarField[(vert_stack[vert_stack_idx][0]+1)*nPointsInSlice + y*nPointsInXDirection + vert_stack[vert_stack_idx][1]];

                var mu2;
                mu2 = (this.m_tIsoLevel - val1)/(val2 - val1);
                var i0_2 = p1x + mu2*(p2x - p1x);
                var i2_2 = p1z + mu2*(p2z - p1z);

                this.chicken_vertices.push(i0_2);
                this.chicken_vertices.push(py);
                this.chicken_vertices.push(i2_2);

                if(this.vert_edge[vert_stack[vert_stack_idx][0]][vert_stack[vert_stack_idx][1]]==1){
                    this.vert_edge[vert_stack[vert_stack_idx][0]][vert_stack[vert_stack_idx][1]]=0;
                    this.ScanAdjacenciesXZ(vert_stack[vert_stack_idx][1],y,vert_stack[vert_stack_idx][0],0,depth);
                    //depth = depth_orig;
                    depth = depth_orig;
                }
                //vert_stack.pop[vert_stack_idx];
                vert_stack_idx--;
            }
        }
        if(horiz_stack.length>0){

            var horiz_stack_idx = horiz_stack.length - 1;
            p1x = x*this.m_fCellLengthX;
            p1z = z*this.m_fCellLengthZ;
            val1 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x];
            if(vert){
                p2x = x*this.m_fCellLengthX;
                p2z = (z+1)*this.m_fCellLengthZ;
                val2 = this.m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + x];
            }else{
                p2x = (x+1)*this.m_fCellLengthX;
                p2z = z*this.m_fCellLengthZ;
                val2 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x+1];
            }

            var mu;
            mu = (this.m_tIsoLevel - val1)/(val2 - val1);
            var i0 = p1x + mu*(p2x - p1x);
            var i2 = p1z + mu*(p2z - p1z);

            while(horiz_stack_idx>=0){
                this.chicken_vertices.push(i0);
                this.chicken_vertices.push(py);
                this.chicken_vertices.push(i2);

                p1x = horiz_stack[horiz_stack_idx][1]*this.m_fCellLengthX;
                p1z = horiz_stack[horiz_stack_idx][0]*this.m_fCellLengthZ;
                p2x = (horiz_stack[horiz_stack_idx][1]+1)*this.m_fCellLengthX;
                p2z = horiz_stack[horiz_stack_idx][0]*this.m_fCellLengthZ;

                val1 = this.m_ptScalarField[horiz_stack[horiz_stack_idx][0]*nPointsInSlice + y*nPointsInXDirection + horiz_stack[horiz_stack_idx][1]];
                val2 = this.m_ptScalarField[horiz_stack[horiz_stack_idx][0]*nPointsInSlice + y*nPointsInXDirection + horiz_stack[horiz_stack_idx][1]+1];

                var mu2;
                mu2 = (this.m_tIsoLevel - val1)/(val2 - val1);
                var i0_2 = p1x + mu2*(p2x - p1x);
                var i2_2 = p1z + mu2*(p2z - p1z);

                this.chicken_vertices.push(i0_2);
                this.chicken_vertices.push(py);
                this.chicken_vertices.push(i2_2);

                if(this.horiz_edge[horiz_stack[horiz_stack_idx][0]][horiz_stack[horiz_stack_idx][1]]==1){
                    this.horiz_edge[horiz_stack[horiz_stack_idx][0]][horiz_stack[horiz_stack_idx][1]]=0;
                    this.ScanAdjacenciesXZ(horiz_stack[horiz_stack_idx][1],y,horiz_stack[horiz_stack_idx][0],1,depth);
                    depth = depth_orig;
                    //depth = depth_orig;
                }
                //horiz_stack.pop[vert_stack_idx];
                horiz_stack_idx--;
            }
        }
    }
}

CIsoSurface.prototype.ScanAdjacencies = function(x, y, z, way, depth){

    var vert_stack  = [];
    var horiz_stack = [];

    var depth_orig = depth;
    if(depth>200){
        //console.log("Warning, returning because of potentially bad deep recursion");
        return;
    }
    depth++;

    var pz = z*this.m_fCellLengthZ;
    var p1x, p1y, p2x,p2y;
    var val1, val2;

    var vert = 0;
    var horiz = 0;
    var drawn = 0;

    var nPointsInXDirection = (this.m_nCellsX + 1);
    var nPointsInSlice = nPointsInXDirection*(this.m_nCellsY + 1);

    var intersection;

    if(way==1){
        horiz = 1;
    } else {
        vert = 1;
    }

    if(x>0&&y>0&&x<this.m_nCellsX-1&&y<this.m_nCellsY-1){

        if(vert){
            if(this.vert_edge[y][x+1]==1)
                vert_stack.push([y,x+1]);
            if(this.vert_edge[y][x-1]==1)
                vert_stack.push([y,x-1]);
            if(this.horiz_edge[y][x]==1)
                horiz_stack.push([y,x]);
            if(this.horiz_edge[y+1][x]==1)
                horiz_stack.push([y+1,x]);
            if(this.horiz_edge[y][x-1]==1)
                horiz_stack.push([y,x-1]);
            if(this.horiz_edge[y+1][x-1]==1)
                horiz_stack.push([y+1,x-1]);
        }else{
            if(this.horiz_edge[y+1][x]==1)
                horiz_stack.push([y+1,x]);
            if(this.horiz_edge[y-1][x]==1)
                horiz_stack.push([y-1,x]);
            if(this.vert_edge[y][x]==1)
                vert_stack.push([y,x]);
            if(this.vert_edge[y-1][x]==1)
                vert_stack.push([y-1,x]);
            if(this.vert_edge[y][x+1]==1)
                vert_stack.push([y,x+1]);
            if(this.vert_edge[y-1][x+1]==1)
                vert_stack.push([y-1,x+1]);
        }

        if(vert_stack.length>0){
            var vert_stack_idx = vert_stack.length - 1;

            p1x = x*this.m_fCellLengthX;
            p1y = y*this.m_fCellLengthY;
            val1 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x];
            if(vert){
                p2x = x*this.m_fCellLengthX;
                p2y = (y+1)*this.m_fCellLengthY;
                val2 = this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x];
            }else{
                p2x = (x+1)*this.m_fCellLengthX;
                p2y = y*this.m_fCellLengthY;
                val2 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x+1];
            }

            var mu;
            mu = (this.m_tIsoLevel - val1)/(val2 - val1);
            var i0 = p1x + mu*(p2x - p1x);
            var i1 = p1y + mu*(p2y - p1y);

            while(vert_stack_idx>=0){
                this.chicken_vertices.push(i0);
                this.chicken_vertices.push(i1);
                this.chicken_vertices.push(pz);

                p1x = vert_stack[vert_stack_idx][1]*this.m_fCellLengthX;
                p1y = vert_stack[vert_stack_idx][0]*this.m_fCellLengthY;
                p2x = vert_stack[vert_stack_idx][1]*this.m_fCellLengthX;
                p2y = (vert_stack[vert_stack_idx][0]+1)*this.m_fCellLengthY;

                val1 = this.m_ptScalarField[z*nPointsInSlice + vert_stack[vert_stack_idx][0]*nPointsInXDirection + vert_stack[vert_stack_idx][1]];
                val2 = this.m_ptScalarField[z*nPointsInSlice + (vert_stack[vert_stack_idx][0]+1)*nPointsInXDirection + vert_stack[vert_stack_idx][1]];

                var mu2;
                mu2 = (this.m_tIsoLevel - val1)/(val2 - val1);
                var i0_2 = p1x + mu2*(p2x - p1x);
                var i1_2 = p1y + mu2*(p2y - p1y);

                this.chicken_vertices.push(i0_2);
                this.chicken_vertices.push(i1_2);
                this.chicken_vertices.push(pz);

                if(this.vert_edge[vert_stack[vert_stack_idx][0]][vert_stack[vert_stack_idx][1]]==1){
                    this.vert_edge[vert_stack[vert_stack_idx][0]][vert_stack[vert_stack_idx][1]]=0;
                    this.ScanAdjacencies(vert_stack[vert_stack_idx][1],vert_stack[vert_stack_idx][0],z,0,depth);
                    //depth = depth_orig;
                    depth = depth_orig;
                }
                //vert_stack.pop();
                vert_stack_idx--;
            }
        }
        if(horiz_stack.length>0){

            var horiz_stack_idx = horiz_stack.length - 1;
            p1x = x*this.m_fCellLengthX;
            p1y = y*this.m_fCellLengthY;
            val1 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x];
            if(vert){
                p2x = x*this.m_fCellLengthX;
                p2y = (y+1)*this.m_fCellLengthY;
                val2 = this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x];
            }else{
                p2x = (x+1)*this.m_fCellLengthX;
                p2y = y*this.m_fCellLengthY;
                val2 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x+1];
            }
            mu = (this.m_tIsoLevel - val1)/(val2 - val1);
            var i0 = p1x + mu*(p2x - p1x);
            var i1 = p1y + mu*(p2y - p1y);

            while(horiz_stack_idx>=0){
                this.chicken_vertices.push(i0);
                this.chicken_vertices.push(i1);
                this.chicken_vertices.push(pz);

                p1x = horiz_stack[horiz_stack_idx][1]*this.m_fCellLengthX;
                p1y = horiz_stack[horiz_stack_idx][0]*this.m_fCellLengthY;
                p2x = (horiz_stack[horiz_stack_idx][1]+1)*this.m_fCellLengthX;
                p2y = horiz_stack[horiz_stack_idx][0]*this.m_fCellLengthY;

                val1 = this.m_ptScalarField[z*nPointsInSlice + horiz_stack[horiz_stack_idx][0]*nPointsInXDirection + horiz_stack[horiz_stack_idx][1]];
                val2 = this.m_ptScalarField[z*nPointsInSlice + horiz_stack[horiz_stack_idx][0]*nPointsInXDirection + horiz_stack[horiz_stack_idx][1]+1];

                //intersection = this.Interpolate(p1x, p1y, pz, p2x, p2y, pz, val1, val2);
                var mu2;
                mu2 = (this.m_tIsoLevel - val1)/(val2 - val1);
                var i0_2 = p1x + mu2*(p2x - p1x);
                var i1_2 = p1y + mu2*(p2y - p1y);

                this.chicken_vertices.push(i0_2);
                this.chicken_vertices.push(i1_2);
                this.chicken_vertices.push(pz);

                if(this.horiz_edge[horiz_stack[horiz_stack_idx][0]][horiz_stack[horiz_stack_idx][1]]==1){
                    this.horiz_edge[horiz_stack[horiz_stack_idx][0]][horiz_stack[horiz_stack_idx][1]]=0;
                    this.ScanAdjacencies(horiz_stack[horiz_stack_idx][1],horiz_stack[horiz_stack_idx][0],z,1,depth);
                    depth = depth_orig;
                    //depth = depth_orig;
                }
                //horiz_stack.pop();
                horiz_stack_idx--;
            }
        }
    }
    //console.log("");
}

CIsoSurface.prototype.GenerateChickenWire = function(ptScalarField, m_tIsoLevel, nCellsX, nCellsY, nCellsZ, fCellLengthX, fCellLengthY, fCellLengthZ){

    this.chicken_vertices = [];

    this.m_nCellsX = nCellsX;
    this.m_nCellsY = nCellsY;
    this.m_nCellsZ = nCellsZ;
    this.m_fCellLengthX = 1.0*fCellLengthX;
    this.m_fCellLengthY = 1.0*fCellLengthY;
    this.m_fCellLengthZ = 1.0*fCellLengthZ;
    this.m_ptScalarField = ptScalarField;

    var nPointsInXDirection = (this.m_nCellsX + 1);
    var nPointsInSlice = nPointsInXDirection*(this.m_nCellsY + 1);

    this.m_tIsoLevel = 1.0*m_tIsoLevel;

    this.vert_edge = [];
    this.horiz_edge = [];
    for (var y = 0; y < this.m_nCellsY; y++){
        this.vert_edge[y] = [];
        this.horiz_edge[y] = [];
    }

    for (var z = 0; z < this.m_nCellsZ; z++){
        for (var y = 0; y < this.m_nCellsY; y++){
            for (var x = 0; x < this.m_nCellsX; x++) {
                this.vert_edge[y][x] = 0;
                this.horiz_edge[y][x] = 0;
                if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x] > this.m_tIsoLevel){
                    if(this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x] <= this.m_tIsoLevel){
                        this.vert_edge[y][x] = 1;
                    }
                    if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + (x+1)] <= this.m_tIsoLevel){
                        this.horiz_edge[y][x] = 1;
                    }
                }
                if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x] < this.m_tIsoLevel){
                    if(this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x] >= this.m_tIsoLevel){
                        this.vert_edge[y][x] = 1;
                    }
                    if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + (x+1)] >= this.m_tIsoLevel){
                        this.horiz_edge[y][x] = 1;
                    }
                }
            }
        }

        for (var y = 0; y < this.m_nCellsY; y++){
            for (var x = 0; x < this.m_nCellsX; x++){
                if(this.horiz_edge[y][x]==1){
                    this.horiz_edge[y][x]=0;
                    var depth = 0;
                    this.ScanAdjacencies(x,y,z,1,depth);
                }
                if(this.vert_edge[y][x]==1){
                    this.vert_edge[y][x]=0;
                    var depth = 0;
                    this.ScanAdjacencies(x,y,z,0,depth);
                }
            }
        }
    }

    this.vert_edge = [];
    this.horiz_edge = [];

    for (var z = 0; z < this.m_nCellsZ; z++){
        this.vert_edge[z] = [];
        this.horiz_edge[z] = [];
    }

    for (var y = 0; y < this.m_nCellsY; y++){
        for (var z = 0; z < this.m_nCellsZ; z++){
            for (var x = 0; x < this.m_nCellsX; x++) {
                this.vert_edge[z][x] = 0;
                this.horiz_edge[z][x] = 0;
                if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x] > this.m_tIsoLevel){
                    if(this.m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + x] <= this.m_tIsoLevel){
                        this.vert_edge[z][x] = 1;
                    }
                    if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + (x+1)] <= this.m_tIsoLevel){
                        this.horiz_edge[z][x] = 1;
                    }
                }
                if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x] < this.m_tIsoLevel){
                    if(this.m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + x] >= this.m_tIsoLevel){
                        this.vert_edge[z][x] = 1;
                    }
                    if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + (x+1)] >= this.m_tIsoLevel){
                        this.horiz_edge[z][x] = 1;
                    }
                }
            }
        }

        for (var z = 0; z < this.m_nCellsZ; z++){
            for (var x = 0; x < this.m_nCellsX; x++){
                if(this.horiz_edge[z][x]==1){
                    this.horiz_edge[z][x]=0;
                    var depth = 0;
                    this.ScanAdjacenciesXZ(x,y,z,1,depth);
                }
                if(this.vert_edge[z][x]==1){
                    this.vert_edge[z][x]=0;
                    var depth = 0;
                    this.ScanAdjacenciesXZ(x,y,z,0,depth);
                }
            }
        }
    }

    this.vert_edge = [];
    this.horiz_edge = [];

    for (var z = 0; z < this.m_nCellsZ; z++){
        this.vert_edge[z] = [];
        this.horiz_edge[z] = [];
    }

    for (var x = 0; x < this.m_nCellsX; x++) {
        for (var z = 0; z < this.m_nCellsZ; z++){
            for (var y = 0; y < this.m_nCellsY; y++){
                this.vert_edge[z][y] = 0;
                this.horiz_edge[z][y] = 0;
                if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x] > this.m_tIsoLevel){
                    if(this.m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + x] <= this.m_tIsoLevel){
                        this.vert_edge[z][y] = 1;
                    }
                    if(this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x] <= this.m_tIsoLevel){
                        this.horiz_edge[z][y] = 1;
                    }
                }
                if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x] < this.m_tIsoLevel){
                    if(this.m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + x] >= this.m_tIsoLevel){
                        this.vert_edge[z][y] = 1;
                    }
                    if(this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x] >= this.m_tIsoLevel){
                        this.horiz_edge[z][y] = 1;
                    }
                }
            }
        }

        for (var z = 0; z < this.m_nCellsZ; z++){
            for (var y = 0; y < this.m_nCellsY; y++){
                if(this.horiz_edge[z][y]==1){
                    this.horiz_edge[z][y]=0;
                    var depth = 0;
                    this.ScanAdjacenciesYZ(x,y,z,1,depth);
                }
                if(this.vert_edge[z][y]==1){
                    this.vert_edge[z][y]=0;
                    var depth = 0;
                    this.ScanAdjacenciesYZ(x,y,z,0,depth);
                }
            }
        }
    }

    this.vert_edge = [];
    this.horiz_edge = [];

}

CIsoSurface.prototype.GenerateContours = function(ptScalarField, min_tIsoLevel, max_tIsoLevel, n_contours, nCellsX, nCellsY, nCellsZ, fCellLengthX, fCellLengthY, fCellLengthZ, origin, plane_in, noTransform){

    // noTransform makes sense only for a 2D renderer. We want to not have to transform again.
    this.chicken_vertices = [];

    var ox = origin[0];
    var oy = origin[1];
    var oz = origin[2];

    var scalarfield;

    var ncellsx, ncellsy;

    var plane =      PlaneNormalized(plane_in);
    var plane_orig = PlaneNormalized(plane_in);

    plane[3] = (plane[3] + plane[0]*ox + plane[1]*oy + plane[2]*oz);

    var intersections = [];
    var p1, p2;
    var lx,ly;

    p1 = [0.0,0.0,0.0];
    p2 = [nCellsX*fCellLengthX,0.0,0.0]; 
    intersections.push(line_plane_intersection(plane,p1,p2));

    p1 = [0.0, nCellsY*fCellLengthY,  0.0]; 
    p2 = [nCellsX*fCellLengthX, nCellsY*fCellLengthY, 0.0]; 
    intersections.push(line_plane_intersection(plane,p1,p2));

    p1 = [0.0, nCellsY*fCellLengthY,  nCellsZ*fCellLengthZ]; 
    p2 = [nCellsX*fCellLengthX,  nCellsY*fCellLengthY,  nCellsZ*fCellLengthZ]; 
    intersections.push(line_plane_intersection(plane,p1,p2));

    p1 = [0.0, 0.0,  nCellsZ*fCellLengthZ]; 
    p2 = [nCellsX*fCellLengthX,  0.0,  nCellsZ*fCellLengthZ]; 
    intersections.push(line_plane_intersection(plane,p1,p2));

    p1 = [0.0,0.0,0.0];
    p2 = [0.0,  0.0,  nCellsZ*fCellLengthZ]; 
    intersections.push(line_plane_intersection(plane,p1,p2));

    p1 = [0.0,  nCellsY*fCellLengthY,  0.0];
    p2 = [0.0,  nCellsY*fCellLengthY,  nCellsZ*fCellLengthZ]; 
    intersections.push(line_plane_intersection(plane,p1,p2));

    p1 = [nCellsX*fCellLengthX,  nCellsY*fCellLengthY,  0.0];
    p2 = [nCellsX*fCellLengthX,  nCellsY*fCellLengthY,  nCellsZ*fCellLengthZ]; 
    intersections.push(line_plane_intersection(plane,p1,p2));

    p1 = [nCellsX*fCellLengthX,  0.0,  0.0];
    p2 = [nCellsX*fCellLengthX,  0.0,  nCellsZ*fCellLengthZ]; 
    intersections.push(line_plane_intersection(plane,p1,p2));

    p1 = [0.0, 0.0, 0.0]; 
    p2 = [0.0, nCellsY*fCellLengthY, 0.0]; 
    intersections.push(line_plane_intersection(plane,p1,p2));

    p1 = [0.0, 0.0, nCellsZ*fCellLengthZ]; 
    p2 = [0.0, nCellsY*fCellLengthY, nCellsZ*fCellLengthZ]; 
    intersections.push(line_plane_intersection(plane,p1,p2));

    p1 = [nCellsX*fCellLengthX, 0.0,  nCellsZ*fCellLengthZ]; 
    p2 = [nCellsX*fCellLengthX, nCellsY*fCellLengthY,  nCellsZ*fCellLengthZ];   
    intersections.push(line_plane_intersection(plane,p1,p2));

    p1 = [nCellsX*fCellLengthX, 0.0,  0.0]; 
    p2 = [nCellsX*fCellLengthX, nCellsY*fCellLengthY, 0.0]; 
    intersections.push(line_plane_intersection(plane,p1,p2));

    //console.log(intersections);

    var good_intersections = [];
    var good_intersections_norm = [];
    for(var i=0;i<12;i++){
        if(Math.abs(intersections[i][0])<1e+7&&Math.abs(intersections[i][1])<1e+7&&Math.abs(intersections[i][2])<1e+7&&Math.abs(intersections[i][3])<1e+7){
            //console.log("Checking "+i);
            //console.log(intersections[i]);
            if(intersections[i][0]>=0.0&&intersections[i][0]<=nCellsX*fCellLengthX){
                //console.log("pass x");
                if(intersections[i][1]>=0.0&&intersections[i][1]<=nCellsY*fCellLengthY){
                    //console.log("pass y");
                    if(intersections[i][2]>=0.0&&intersections[i][2]<=nCellsZ*fCellLengthZ){
                        //console.log("pass z");
                        var norm_intersection = VectorNormalized(intersections[i]);
                        if(good_intersections_norm.length<1){
                            //console.log("Adding "+i); 
                            good_intersections.push(intersections[i]);
                            good_intersections_norm.push(norm_intersection);
                        }else{
                            var had_before = false;
                            for(jj=0;jj<good_intersections_norm.length;jj++){
                                if(Math.abs(DotProduct(good_intersections_norm[jj],norm_intersection))>.9){
                                    if(VectorLength(intersections[i])>VectorLength(good_intersections_norm[jj])){
                                        good_intersections[jj] = intersections[i];
                                        good_intersections_norm[jj] = norm_intersection;
                                    }
                                    had_before = true;
                                }
                            }
                            if(!had_before){
                                //console.log("Adding "+i); 
                                good_intersections.push(intersections[i]);
                                good_intersections_norm.push(norm_intersection);
                            }
                        }
                    }
                }
            }
        }
    }

    if(good_intersections.length<4){
        console.log("Doesn't intersect four lines, returning\n");
        return;
    }

    if(good_intersections.length>4){
        while(good_intersections.length>4){
            // Find a non redundant set ...
            var max_overlap = 0.0;
            var max_overlap_ind = -1;
            for(ii=0;ii<good_intersections.length;ii++){
                var this_norm_ii = VectorNormalized(good_intersections[ii]);
                for(jj=0;jj<good_intersections.length;jj++){
                    if(ii!=jj){
                        var this_norm_jj = VectorNormalized(good_intersections[jj]);
                        var overlap = Math.abs(DotProduct(this_norm_ii,this_norm_jj));
                        if(overlap>max_overlap){
                            max_overlap = overlap;
                            max_overlap_ind = jj;
                        }
                    }
                }
            }
            good_intersections.splice(max_overlap_ind,1);
        }
    }

    var PIBY2 = Math.PI * 2;
    var xprime;
    var yprime;
    var xprime_orig;
    var yprime_orig;
    var diff_from_origin = [0,0,0,1];
    if(good_intersections.length>=4){
        var angle = Angle(good_intersections[0],good_intersections[1],good_intersections[2]);
        angle += Angle(good_intersections[1],good_intersections[2],good_intersections[3]);
        angle += Angle(good_intersections[2],good_intersections[3],good_intersections[0]);
        angle += Angle(good_intersections[3],good_intersections[0],good_intersections[1]);
        //std::cout << "total internal angles of quad: " << angle << std::endl;
        if(Math.abs(PIBY2-angle)>0.0001){ // Could be crossed, swap vertices.
            //std::cout << "Crossing!!!!!!!!!!!\n";
            intersections = [];
            intersections.push(good_intersections[3]);
            intersections.push(good_intersections[1]);
            intersections.push(good_intersections[0]);
            intersections.push(good_intersections[2]);
            good_intersections = intersections;
            angle = Angle(good_intersections[0],good_intersections[1],good_intersections[2]);
            angle += Angle(good_intersections[1],good_intersections[2],good_intersections[3]);
            angle += Angle(good_intersections[2],good_intersections[3],good_intersections[0]);
            angle += Angle(good_intersections[3],good_intersections[0],good_intersections[1]);
            //std::cout << "After perm total internal angles of quad: " << angle << std::endl;
            if(Math.abs(PIBY2-angle)>0.0001){ // Bigger problem.
                intersections = [];
                intersections.push(good_intersections[0]);
                intersections.push(good_intersections[1]);
                intersections.push(good_intersections[3]);
                intersections.push(good_intersections[2]);
                good_intersections = intersections;
                angle = Angle(good_intersections[0],good_intersections[1],good_intersections[2]);
                angle += Angle(good_intersections[1],good_intersections[2],good_intersections[3]);
                angle += Angle(good_intersections[2],good_intersections[3],good_intersections[0]);
                angle += Angle(good_intersections[3],good_intersections[0],good_intersections[1]);
                //std::cout << "After perm(2) total internal angles of quad: " << angle << std::endl;
                if(Math.abs(PIBY2-angle)>0.0001) // Bigger problem.
                    return;
            }
        }

        var plane_mid_point = VectorArrayMidPoint(good_intersections);

        diff_from_origin  = VectorDiff(plane_mid_point , [-ox,-oy,-oz]);

        var l1 = VectorLength(VectorDiff(good_intersections[0] , good_intersections[1]));
        var l2 = VectorLength(VectorDiff(good_intersections[1] , good_intersections[2]));
        var l3 = VectorLength(VectorDiff(good_intersections[2] , good_intersections[3]));
        var l4 = VectorLength(VectorDiff(good_intersections[3] , good_intersections[0]));
        //std::cout << "Length of sides: " << l1 << " " << l2 << " " << l3 << " " << l4 << std::endl;
        if(l1<(nCellsX*fCellLengthX)/10||l2<(nCellsX*fCellLengthX)/10||l3<(nCellsX*fCellLengthX)/10||l4<(nCellsX*fCellLengthX)/10){
            console.log("New grid too small: " ,(nCellsX*fCellLengthX)/10);
            return;
        }
        lx = l1;
        ly = l2;
        if(l3<lx)
            lx = l3;
        if(l4<ly)
            ly = l4;

        xprime = VectorDiff(good_intersections[3] , good_intersections[0]);
        yprime = VectorDiff(good_intersections[1] , good_intersections[0]);

        xprime_orig = VectorCopy(xprime);
        xprime = VectorNormalized(xprime);

        var plane_norm = VectorNormalized(PlaneGetNormal(plane));

        var new_y_prime = CrossProduct(xprime,plane_norm);
        if(DotProduct(new_y_prime,yprime)<0.0) new_y_prime = VectorNegate(new_y_prime);

        good_intersections[1] = VectorSum(good_intersections[0] , FloatVectorMult(VectorLength(yprime),new_y_prime));
        yprime = VectorDiff(good_intersections[1],good_intersections[0]);

        yprime_orig = VectorCopy(yprime);

        yprime = VectorNormalized(yprime);

        var ncellsx_f = lx/fCellLengthX;
        var ncellsy_f = ly/fCellLengthY;
        ncellsx = parseInt(Math.floor(ncellsx_f+0.5));
        ncellsy = parseInt(Math.floor(ncellsy_f+0.5));

        scalarfield = [];
        for (var y = 0; y < ncellsy; y++){
            var yfrac = 1.0*y/ncellsy;
            var py =  VectorSum(FloatVectorMult(yfrac,good_intersections[1]),FloatVectorMult((1.0-yfrac),good_intersections[0]));
            var pyp = VectorSum(FloatVectorMult(yfrac,good_intersections[2]),FloatVectorMult((1.0-yfrac),good_intersections[3]));
            for (var x = 0; x < ncellsx; x++) {
                var xfrac = 1.0*x/ncellsx;
                var p = VectorSum(py,FloatVectorMult(xfrac,VectorDiff(pyp,py)));
                var fx = parseInt(p[0]/fCellLengthX);
                var fy = parseInt(p[1]/fCellLengthY);
                var fz = parseInt(p[2]/fCellLengthZ);

                var intp = Math.floor(p[0]/fCellLengthX);
                var xinterp = p[0]/fCellLengthX-intp;
                var intp = Math.floor(p[1]/fCellLengthY);
                var yinterp = p[1]/fCellLengthY-intp;
                var intp = Math.floor(p[2]/fCellLengthZ);
                var zinterp = p[2]/fCellLengthZ-intp;

                var ind1 = fz*(nCellsX+1)*(nCellsY+1) + fy*(nCellsX+1) + fx;
                var ind2 = fz*(nCellsX+1)*(nCellsY+1) + fy*(nCellsX+1) + fx+1;
                var ind3 = fz*(nCellsX+1)*(nCellsY+1) + (fy+1)*(nCellsX+1) + fx;
                var ind4 = fz*(nCellsX+1)*(nCellsY+1) + (fy+1)*(nCellsX+1) + fx+1;
                var ind5 = (fz+1)*(nCellsX+1)*(nCellsY+1) + fy*(nCellsX+1) + fx;
                var ind6 = (fz+1)*(nCellsX+1)*(nCellsY+1) + fy*(nCellsX+1) + fx+1;
                var ind7 = (fz+1)*(nCellsX+1)*(nCellsY+1) + (fy+1)*(nCellsX+1) + fx;
                var ind8 = (fz+1)*(nCellsX+1)*(nCellsY+1) + (fy+1)*(nCellsX+1) + fx+1;
                if(ind1>=(nCellsX+1)*(nCellsY+1)*(nCellsZ+1)) continue;
                if(ind2>=(nCellsX+1)*(nCellsY+1)*(nCellsZ+1)) continue;
                if(ind3>=(nCellsX+1)*(nCellsY+1)*(nCellsZ+1)) continue;
                if(ind4>=(nCellsX+1)*(nCellsY+1)*(nCellsZ+1)) continue;
                if(ind5>=(nCellsX+1)*(nCellsY+1)*(nCellsZ+1)) continue;
                if(ind6>=(nCellsX+1)*(nCellsY+1)*(nCellsZ+1)) continue;
                if(ind7>=(nCellsX+1)*(nCellsY+1)*(nCellsZ+1)) continue;
                if(ind8>=(nCellsX+1)*(nCellsY+1)*(nCellsZ+1)) continue;
                if(y*ncellsx+x>=(ncellsx)*(ncellsy)) continue;
                var val1 = ptScalarField[ind1];
                var val2 = ptScalarField[ind2];
                var valx1 = val1+ xinterp*(val2-val1);
                val1 = ptScalarField[ind3];
                val2 = ptScalarField[ind4];
                var valx2 = val1+ xinterp*(val2-val1);
                var valxy1 = valx1 + yinterp*(valx2-valx1);
                val1 = ptScalarField[ind5];
                val2 = ptScalarField[ind6];
                valx1 = val1+ xinterp*(val2-val1);
                val1 = ptScalarField[ind7];
                val2 = ptScalarField[ind8];
                valx2 = val1+ xinterp*(val2-val1);
                var valxy2 = valx1 + yinterp*(valx2-valx1);
                var valxyz = valxy1 + zinterp*(valxy2-valxy1);
                scalarfield[y*ncellsx+x] = valxyz;
            }
        }
    }else{
        console.log("Doesn't intersect four lines, returning");
        return;
    }

    var z = 0; 

    this.m_nCellsX = nCellsX;
    this.m_nCellsY = nCellsY;
    this.m_nCellsZ = nCellsZ;
    this.m_fCellLengthX = fCellLengthX;
    this.m_fCellLengthY = fCellLengthY;
    this.m_fCellLengthZ = fCellLengthZ;
    this.m_ptScalarField = ptScalarField;

    this.m_nCellsX = ncellsx-1;
    this.m_nCellsY = ncellsy-1;
    this.m_nCellsZ = 1;
    this.m_ptScalarField = scalarfield;

    var nPointsInXDirection = (this.m_nCellsX + 1);
    var nPointsInSlice = nPointsInXDirection*(this.m_nCellsY + 1);

    this.vert_edge = [];
    this.horiz_edge = [];
    for (var y = 0; y < this.m_nCellsY; y++){
        this.vert_edge[y] = [];
        this.horiz_edge[y] = [];
    }

    this.m_tIsoLevel = min_tIsoLevel;

    //int n_contours = (int)((max_tIsoLevel-min_tIsoLevel)/contour_sub_divide);
    var contour_sub_divide = (max_tIsoLevel-min_tIsoLevel)/n_contours;

    //console.log(contour_sub_divide);
    //console.log(ptScalarField);
    //console.log(scalarfield);

    for(var icont=0;icont<n_contours;icont++){

        this.m_tIsoLevel = min_tIsoLevel + icont* contour_sub_divide;
        for (var y = 0; y < this.m_nCellsY; y++){
            for (var x = 0; x < this.m_nCellsX; x++) {
                this.vert_edge[y][x] = 0;
                this.horiz_edge[y][x] = 0;
                if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x] > this.m_tIsoLevel){
                    if(this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x] <= this.m_tIsoLevel){
                        this.vert_edge[y][x] = 1;
                    }
                    if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + (x+1)] <= this.m_tIsoLevel){
                        this.horiz_edge[y][x] = 1;
                    }
                }
                if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x] < this.m_tIsoLevel){
                    if(this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x] >= this.m_tIsoLevel){
                        this.vert_edge[y][x] = 1;
                    }
                    if(this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + (x+1)] >= this.m_tIsoLevel){
                        this.horiz_edge[y][x] = 1;
                    }
                }
            }
        }

        for (var y = 0; y < this.m_nCellsY; y++){
            for (var x = 0; x < this.m_nCellsX; x++){
                if(this.horiz_edge[y][x]==1){
                    var p1x = x*this.m_fCellLengthX;
                    var p1y = y*this.m_fCellLengthY;
                    var pz = z*this.m_fCellLengthZ;
                    var val1 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x];
                    var p2x = (x+1)*this.m_fCellLengthX;
                    var p2y = y*this.m_fCellLengthY;
                    var val2 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x+1];
                    var intersection = this.Interpolate(p1x, p1y, pz, p2x, p2y, pz, val1, val2);
                    var depth = 0;
                    this.ScanAdjacencies(x,y,z,1,depth);
                    this.horiz_edge[y][x]=0;
                }
                if(this.vert_edge[y][x]==1){
                    //std::cout << "Vert: " << x << " " << y << std::endl;
                    var p1x = x*this.m_fCellLengthX;
                    var p1y = y*this.m_fCellLengthY;
                    var pz = z*this.m_fCellLengthZ;
                    var val1 = this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x];
                    var p2x = x*this.m_fCellLengthX;
                    var p2y = (y+1)*this.m_fCellLengthY;
                    var val2 = this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x];
                    var intersection = this.Interpolate(p1x, p1y, pz, p2x, p2y, pz, val1, val2);
                    var depth = 0;
                    this.ScanAdjacencies(x,y,z,0,depth);
                    this.vert_edge[y][x]=0;
                }
            }
        }
    }

    for (var y = 0; y < this.m_nCellsY; y++){
        this.vert_edge[y] = []; 
        this.horiz_edge[y] = [];
    }
    this.vert_edge = [];
    this.horiz_edge = [];

    // At this point we should transform the coords to the right plane ... 

    if(noTransform){
        return;
    }

    //std::cout << "Original Plane: " << plane_orig.get_A() << " " << plane_orig.get_B() << " " << plane_orig.get_C() << " " << plane_orig.get_D() << "\n";
    plane = PlaneNormalized(plane);
    //std::cout << "Normalized/moved Plane: " << plane.get_A() << " " << plane.get_B() << " " << plane.get_C() << " " << plane.get_D() << "\n";
    //std::cout << "Plane Centre: " << lx/2.0 << " " << ly/2.0 << std::endl;
    //std::cout << "Origin: " << ox << " " << oy << " " << oz << std::endl;

    var plane_norm = PlaneGetNormal(plane);

    var xaxis = [1.0, 0.0, 0.0, 1.0];
    var yaxis = [0.0, 1.0, 0.0, 1.0];
    var zaxis = [0.0, 0.0, 1.0, 1.0];
    var cross = CrossProduct(plane_norm,zaxis);

    if(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]>1.0e-6){

        //double sign = plane.get_C();
        var length = VectorLength(plane_norm);
        var angle = Math.acos(Math.abs(plane[2]/length));

        if(plane[2]>0.0) angle = -angle;

        //std::cout << "Will rotate by: " << angle << " radians about (" << cross[0] << ", " << cross[1] << ", " << cross[2] << ")\n";
        var rotquat = Quat(cross[0],cross[1],cross[2],1,-angle*180.0/Math.PI);

        var invQuat = quat4.create();
        quat4.inverse(rotquat,invQuat);
        var rotmat =  quat4.toMat4(invQuat);

        var pos = [];
        pos[3] = 1.0;

        var result = vec3.create();
        var x_scale = (1.0*ncellsx)/(1.0*nCellsX);
        var y_scale = (1.0*ncellsy)/(1.0*nCellsY);

        var fRotMat = quat4.toMat4(rotquat);
        var xprime_rot = vec3.create();
        var yprime_rot = vec3.create();

        var xprime_v = vec3.create(xprime);
        var yprime_v = vec3.create(yprime);

        xprime_rot = RotMxV(fRotMat,xprime_v);
        yprime_rot = RotMxV(fRotMat,yprime_v);

        //std::cout << "x: " << xprime_rot << "\n"; std::cout.flush();
        //std::cout << "y: " << yprime_rot << "\n"; std::cout.flush();

        //std::cout << Cartesian::DotProduct(xprime_rot,xaxis) << " " << Cartesian::DotProduct(yprime_rot,yaxis) << "\n";
        //std::cout << " " << Cartesian::DotProduct(xprime_rot,yprime_rot) << " " << Cartesian::DotProduct(xaxis,yaxis) << "\n";
        if(Math.abs(DotProduct(xprime_rot,xaxis))<0.8&&Math.abs(DotProduct(yprime_rot,xaxis))>0.8) {
            //std::cout << "Swapping x and y\n"; std::cout.flush();
            var ctmp = VectorCopy(xprime);
            xprime = VectorCopy(yprime);
            yprime = VectorCopy(ctmp);
            var z_ang = 180.0/M_PI * Math.acos(DotProduct(xprime_rot,yaxis));
            //std::cout << z_ang << "\n";
            var rotz = Quat(0,0,1,1,-90+z_ang);
            quat4.multiply(rotquat,rotz);
            quat4.inverse(rotquat,invQuat);
            rotmat =  quat4.toMat4(invQuat);
        } else if(Math.abs(DotProduct(xprime_rot,xaxis))-1.0<1e-7){
            var rotz = Quat(0,0,1,1,0.0);
            quat4.multiply(rotquat,rotz);
            quat4.inverse(rotquat,invQuat);
            rotmat =  quat4.toMat4(invQuat);
            //std::cout << "rotating by " << 0.0 << " about z\n";
        } else {
            var z_ang = 180.0/M_PI * Math.acos(DotProduct(xprime_rot,xaxis));
            if(z_ang>90) z_ang = 180 - z_ang;
            var rotz = Quat(0,0,1,1,z_ang);
            quat4.multiply(rotquat,rotz);
            quat4.inverse(rotquat,invQuat);
            rotmat =  quat4.toMat4(invQuat);
            //std::cout << "rotating by " << z_ang << " about z\n";
        }

        if(DotProduct(yprime_rot,yaxis)<0.0) x_scale = -x_scale;
        if(DotProduct(xprime_rot,xaxis)<0.0) y_scale = -y_scale;

        //std::cout << "x_scale: " << x_scale << "\n"; std::cout.flush();
        //std::cout << "y_scale: " << y_scale << "\n"; std::cout.flush();

        for(var i=0;i<this.chicken_vertices.length;i+=3){
            //std::cout << "Before adding origin: " << chicken_vertices[i] << " " << chicken_vertices[i+1] << " " << chicken_vertices[i+2] << std::endl;
            this.chicken_vertices[i]   -= lx/2.0;
            this.chicken_vertices[i+1] -= ly/2.0;
            //std::cout << "After adding origin: " << chicken_vertices[i] << " " << chicken_vertices[i+1] << " " << chicken_vertices[i+2] << std::endl;
            pos[0] = this.chicken_vertices[i]*(this.m_fCellLengthX*nCellsX)/lx * y_scale; // Note swap of x and y!
            pos[1] = this.chicken_vertices[i+1]*(this.m_fCellLengthY*nCellsY)/ly * x_scale; // SCALINGS MAY BE DODGY !!!!!
            pos[2] = this.chicken_vertices[i+2];

            var pos_v = vec3.create(pos);
            result = RotMxV(rotmat,pos_v);

            //std::cout << "After rotating: " << result[0] << " " << result[1] << " " << result[2] << std::endl;
            this.chicken_vertices[i]   = result[0]+diff_from_origin[0];
            this.chicken_vertices[i+1] = result[1]+diff_from_origin[1];
            this.chicken_vertices[i+2] = result[2]+diff_from_origin[2];
            //std::cout << "After adding plane: " << chicken_vertices[i] << " " << chicken_vertices[i+1] << " " << chicken_vertices[i+2] << std::endl;
        }
    }else{
        for(var i=0;i<this.chicken_vertices.length;i+=3){
            this.chicken_vertices[i]   += ox;
            this.chicken_vertices[i+1] += oy;
            this.chicken_vertices[i+2] -= plane_orig[3]*plane_orig[2];
        }
    }
    //console.log(this.chicken_vertices);
}

CIsoSurface.prototype.GenerateContoursFromDiv = function(divid, tIsoLevel, nconts, nCellsX, nCellsY, nCellsZ, fCellLengthX, fCellLengthY, fCellLengthZ, origin, plane, noTransform){
    var theIsland = document.getElementById(divid).firstChild;
    var ptScalarField = JSON.parse(theIsland.data);
    this.GenerateContours(ptScalarField, 0.0, tIsoLevel, nconts, nCellsX, nCellsY, nCellsZ, fCellLengthX, fCellLengthY, fCellLengthZ, origin, plane, noTransform);

}

CIsoSurface.prototype.GenerateSurfaceFromDiv = function(divid, tIsoLevel, nCellsX, nCellsY, nCellsZ, fCellLengthX, fCellLengthY, fCellLengthZ)
{
    //console.log("Let's load the data ...");
    var theIsland = document.getElementById(divid).firstChild;
    //console.log(JSON.parse(theIsland.data));
    //console.log(JSON.parse(theIsland.data).length);
    var ptScalarField = JSON.parse(theIsland.data);
    this.GenerateSurface(ptScalarField, tIsoLevel, nCellsX, nCellsY, nCellsZ, fCellLengthX, fCellLengthY, fCellLengthZ);
}

CIsoSurface.prototype.GenerateSurfacePartial = function(ptScalarField, tIsoLevel, nCellsX, nCellsY, nCellsZ, fCellLengthX, fCellLengthY, fCellLengthZ, baseCellX, baseCellY, baseCellZ, upperCellX, upperCellY, upperCellZ, doChickenWire)
{
    var m_baseCellX = 0;
    var m_baseCellY = 0;
    var m_baseCellZ = 0;
    var m_upperCellX = nCellsX;
    var m_upperCellY = nCellsY;
    var m_upperCellZ = nCellsZ;

    var nPointsInXDirection = (nCellsX + 1);
    var nPointsInSlice = nPointsInXDirection*(nCellsY + 1);

    if(typeof upperCellX!=="undefined"&&typeof upperCellY!=="undefined"&&typeof upperCellZ!=="undefined"&&typeof baseCellX!=="undefined"&&typeof baseCellY!=="undefined"&&typeof baseCellZ!=="undefined"){
        //console.log("changing bounds");
        m_baseCellX = baseCellX;
        m_baseCellY = baseCellY;
        m_baseCellZ = baseCellZ;
        m_upperCellX = upperCellX;
        m_upperCellY = upperCellY;
        m_upperCellZ = upperCellZ;
    }

    var new_nPointsInXDirection = (m_upperCellX-m_baseCellX + 1);
    var new_nPointsInSlice = new_nPointsInXDirection*(m_upperCellY-m_baseCellY + 1);

    var newBuffer = [];

    var zNew = 0;
    for (var z = m_baseCellZ; z <= upperCellZ; z++,zNew++){
        var yNew = 0;
        for (var y = m_baseCellY; y <= upperCellY; y++,yNew++){
            var xNew = 0;
            for (var x = m_baseCellX; x <= upperCellX; x++,xNew++) {
                newBuffer[zNew*new_nPointsInSlice + yNew*new_nPointsInXDirection + xNew] = ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x];
            }
        }
    }
    /*
    console.log("The new stuff .......");
    console.log(newBuffer.length);
    console.log(tIsoLevel);
    console.log(m_upperCellX-m_baseCellX);
    console.log(m_upperCellY-m_baseCellY);
    console.log(m_upperCellZ-m_baseCellZ);
    console.log(fCellLengthX);
    console.log(fCellLengthY);
    console.log(fCellLengthZ);
    */
    if(typeof(doChickenWire)!=="undefined" && doChickenWire){
        this.GenerateChickenWire(newBuffer,tIsoLevel, m_upperCellX-m_baseCellX, m_upperCellY-m_baseCellY, m_upperCellZ-m_baseCellZ, fCellLengthX, fCellLengthY, fCellLengthZ);
    }else{
        this.GenerateSurface(newBuffer,tIsoLevel, m_upperCellX-m_baseCellX, m_upperCellY-m_baseCellY, m_upperCellZ-m_baseCellZ, fCellLengthX, fCellLengthY, fCellLengthZ);
    }
}

CIsoSurface.prototype.GenerateSurface = function(ptScalarField, tIsoLevel, nCellsX, nCellsY, nCellsZ, fCellLengthX, fCellLengthY, fCellLengthZ)
{
    if (this.m_bValidSurface)
        this.DeleteSurface();

    var start = new Date().getTime();
    this.m_tIsoLevel = tIsoLevel;
    this.m_nCellsX = nCellsX;
    this.m_nCellsY = nCellsY;
    this.m_nCellsZ = nCellsZ;
    if(this.cell){
        // We'll do all the scaling, etc. in returnVertices.
        this.m_fCellLengthX = 1.0;
        this.m_fCellLengthY = 1.0;
        this.m_fCellLengthZ = 1.0;
    } else {
        this.m_fCellLengthX = fCellLengthX;
        this.m_fCellLengthY = fCellLengthY;
        this.m_fCellLengthZ = fCellLengthZ;
    }
    this.m_ptScalarField = ptScalarField;

    var nPointsInXDirection = (this.m_nCellsX + 1);
    var nPointsInSlice = nPointsInXDirection*(this.m_nCellsY + 1);

    // Generate isosurface.
    for (var z = 0; z < this.m_nCellsZ; z++){
        for (var y = 0; y < this.m_nCellsY; y++){
            for (var x = 0; x < this.m_nCellsX; x++) {
                // Calculate table lookup index from those
                // vertices which are below the isolevel.
                var tableIndex = 0;
                if (this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + x] < this.m_tIsoLevel)
                    tableIndex |= 1;
                if (this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + x] < this.m_tIsoLevel)
                    tableIndex |= 2;
                if (this.m_ptScalarField[z*nPointsInSlice + (y+1)*nPointsInXDirection + (x+1)] < this.m_tIsoLevel)
                    tableIndex |= 4;
                if (this.m_ptScalarField[z*nPointsInSlice + y*nPointsInXDirection + (x+1)] < this.m_tIsoLevel)
                    tableIndex |= 8;
                if (this.m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + x] < this.m_tIsoLevel)
                    tableIndex |= 16;
                if (this.m_ptScalarField[(z+1)*nPointsInSlice + (y+1)*nPointsInXDirection + x] < this.m_tIsoLevel)
                    tableIndex |= 32;
                if (this.m_ptScalarField[(z+1)*nPointsInSlice + (y+1)*nPointsInXDirection + (x+1)] < this.m_tIsoLevel)
                    tableIndex |= 64;
                if (this.m_ptScalarField[(z+1)*nPointsInSlice + y*nPointsInXDirection + (x+1)] < this.m_tIsoLevel)
                    tableIndex |= 128;

                // Now create a triangulation of the isosurface in this
                // cell.
                if (m_edgeTable[tableIndex] != 0) {
                    if (m_edgeTable[tableIndex] & 8) {
                        var pt = this.CalculateIntersection(x, y, z, 3);
                        var id = this.GetEdgeID(x, y, z, 3);
                        this.m_i2pt3idVertices[id] = pt;
                    }
                    if (m_edgeTable[tableIndex] & 1) {
                        var pt = this.CalculateIntersection(x, y, z, 0);
                        var id = this.GetEdgeID(x, y, z, 0);
                        this.m_i2pt3idVertices[id] = pt;
                    }
                    if (m_edgeTable[tableIndex] & 256) {
                        var pt = this.CalculateIntersection(x, y, z, 8);
                        var id = this.GetEdgeID(x, y, z, 8);
                        this.m_i2pt3idVertices[id] = pt;
                    }

                    if (x == this.m_nCellsX - 1) {
                        if (m_edgeTable[tableIndex] & 4) {
                            var pt = this.CalculateIntersection(x, y, z, 2);
                            var id = this.GetEdgeID(x, y, z, 2);
                            this.m_i2pt3idVertices[id] = pt;
                        }
                        if (m_edgeTable[tableIndex] & 2048) {
                            var pt = this.CalculateIntersection(x, y, z, 11);
                            var id = this.GetEdgeID(x, y, z, 11);
                            this.m_i2pt3idVertices[id] = pt;
                        }
                    }
                    if (y == this.m_nCellsY - 1) {
                        if (m_edgeTable[tableIndex] & 2) {
                            var pt = this.CalculateIntersection(x, y, z, 1);
                            var id = this.GetEdgeID(x, y, z, 1);
                            this.m_i2pt3idVertices[id] = pt;
                        }
                        if (m_edgeTable[tableIndex] & 512) {
                            var pt = this.CalculateIntersection(x, y, z, 9);
                            var id = this.GetEdgeID(x, y, z, 9);
                            this.m_i2pt3idVertices[id] = pt;
                        }
                    }
                    if (z == this.m_nCellsZ - 1) {
                        if (m_edgeTable[tableIndex] & 16) {
                            var pt = this.CalculateIntersection(x, y, z, 4);
                            var id = this.GetEdgeID(x, y, z, 4);
                            this.m_i2pt3idVertices[id] = pt;
                        }
                        if (m_edgeTable[tableIndex] & 128) {
                            var pt = this.CalculateIntersection(x, y, z, 7);
                            var id = this.GetEdgeID(x, y, z, 7);
                            this.m_i2pt3idVertices[id] = pt;
                        }
                    }
                    if ((x==this.m_nCellsX - 1) && (y==this.m_nCellsY - 1))
                        if (m_edgeTable[tableIndex] & 1024) {
                            var pt = this.CalculateIntersection(x, y, z, 10);
                            var id = this.GetEdgeID(x, y, z, 10);
                            this.m_i2pt3idVertices[id] = pt;
                        }
                    if ((x==this.m_nCellsX - 1) && (z==this.m_nCellsZ - 1))
                        if (m_edgeTable[tableIndex] & 64) {
                            var pt = this.CalculateIntersection(x, y, z, 6);
                            var id = this.GetEdgeID(x, y, z, 6);
                            this.m_i2pt3idVertices[id] = pt;
                        }
                    if ((y==this.m_nCellsY - 1) && (z==this.m_nCellsZ - 1))
                        if (m_edgeTable[tableIndex] & 32) {
                            var pt = this.CalculateIntersection(x, y, z, 5);
                            var id = this.GetEdgeID(x, y, z, 5);
                            this.m_i2pt3idVertices[id] = pt;
                        }

                    for (var i = 0; m_triTable[tableIndex][i] != -1; i += 3) {
                        var triangle = [-1,-1,-1,-1];
                        var pointID0, pointID1, pointID2;
                        pointID0 = this.GetEdgeID(x, y, z, m_triTable[tableIndex][i]);
                        pointID1 = this.GetEdgeID(x, y, z, m_triTable[tableIndex][i+1]);
                        pointID2 = this.GetEdgeID(x, y, z, m_triTable[tableIndex][i+2]);
                        triangle[0] = pointID0;
                        triangle[1] = pointID1;
                        triangle[2] = pointID2;
                        this.m_trivecTriangles.push(triangle);
                    }
                }
            }
        }
    }

    var end = new Date().getTime();
    var time = end - start;
    console.log("Time to do bulk of recalc maps: "+time*0.001);

    this.RenameVerticesAndTriangles();
    end = new Date().getTime();
    var time = end - start;
    console.log("Time to do RenameVerticesAndTriangles: "+time*0.001);
    //CalculateNormals();
    this.m_bValidSurface = true;
}


 

CIsoSurface.prototype.IsSurfaceValid = function()
{
    return this.m_bValidSurface;
}

CIsoSurface.prototype.DeleteSurface = function(){
    this.m_fCellLengthX = 0;
    this.m_fCellLengthY = 0;
    this.m_fCellLengthZ = 0;
    this.m_nCellsX = 0;
    this.m_nCellsY = 0;
    this.m_nCellsZ = 0;
    this.m_nTriangles = 0;
    this.m_nTriangles_clipped = 0;
    this.m_nNormals = 0;
    this.m_nVertices = 0;
    this.m_ppt3dVertices = [];
    this.m_piTriangleIndices = [];
    this.m_pvec3dNormals = [];
    this.m_ptScalarField = [];
    this.m_tIsoLevel = 0;
    this.m_bValidSurface = false;
    this.m_i2pt3idVertices = {};
    this.m_trivecTriangles = [];
    this.chicken_vertices = [];
}

CIsoSurface.prototype.GetEdgeID = function(nX, nY, nZ, nEdgeNo)
{
	switch (nEdgeNo) {
	case 0:
		return this.GetVertexID(nX, nY, nZ) + 1;
	case 1:
		return this.GetVertexID(nX, nY + 1, nZ);
	case 2:
		return this.GetVertexID(nX + 1, nY, nZ) + 1;
	case 3:
		return this.GetVertexID(nX, nY, nZ);
	case 4:
		return this.GetVertexID(nX, nY, nZ + 1) + 1;
	case 5:
		return this.GetVertexID(nX, nY + 1, nZ + 1);
	case 6:
		return this.GetVertexID(nX + 1, nY, nZ + 1) + 1;
	case 7:
		return this.GetVertexID(nX, nY, nZ + 1);
	case 8:
		return this.GetVertexID(nX, nY, nZ) + 2;
	case 9:
		return this.GetVertexID(nX, nY + 1, nZ) + 2;
	case 10:
		return this.GetVertexID(nX + 1, nY + 1, nZ) + 2;
	case 11:
                return this.GetVertexID(nX + 1, nY, nZ) + 2;
	default:
		// Invalid edge no.
		return -1;
	}
}

CIsoSurface.prototype.GetVertexID = function(nX, nY, nZ)
{
	return 3*(nZ*(this.m_nCellsY + 1)*(this.m_nCellsX + 1) + nY*(this.m_nCellsX + 1) + nX);
}

CIsoSurface.prototype.CalculateIntersection = function(nX, nY, nZ, nEdgeNo)
{
    var x1, y1, z1, x2, y2, z2;
    var v1x = nX, v1y = nY, v1z = nZ;
    var v2x = nX, v2y = nY, v2z = nZ;

    switch (nEdgeNo)
    {
        case 0:
            v2y += 1;
            break;
        case 1:
            v1y += 1;
            v2x += 1;
            v2y += 1;
            break;
        case 2:
            v1x += 1;
            v1y += 1;
            v2x += 1;
            break;
        case 3:
            v1x += 1;
            break;
        case 4:
            v1z += 1;
            v2y += 1;
            v2z += 1;
            break;
        case 5:
            v1y += 1;
            v1z += 1;
            v2x += 1;
            v2y += 1;
            v2z += 1;
            break;
        case 6:
            v1x += 1;
            v1y += 1;
            v1z += 1;
            v2x += 1;
            v2z += 1;
            break;
        case 7:
            v1x += 1;
            v1z += 1;
            v2z += 1;
            break;
        case 8:
            v2z += 1;
            break;
        case 9:
            v1y += 1;
            v2y += 1;
            v2z += 1;
            break;
        case 10:
            v1x += 1;
            v1y += 1;
            v2x += 1;
            v2y += 1;
            v2z += 1;
            break;
        case 11:
            v1x += 1;
            v2x += 1;
            v2z += 1;
            break;
    }

    x1 = v1x*this.m_fCellLengthX;
    y1 = v1y*this.m_fCellLengthY;
    z1 = v1z*this.m_fCellLengthZ;
    x2 = v2x*this.m_fCellLengthX;
    y2 = v2y*this.m_fCellLengthY;
    z2 = v2z*this.m_fCellLengthZ;

    var nPointsInXDirection = (this.m_nCellsX + 1);
    var nPointsInSlice = nPointsInXDirection*(this.m_nCellsY + 1);
    var val1 = this.m_ptScalarField[v1z*nPointsInSlice + v1y*nPointsInXDirection + v1x];
    var val2 = this.m_ptScalarField[v2z*nPointsInSlice + v2y*nPointsInXDirection + v2x];
    var intersection = this.Interpolate(x1, y1, z1, x2, y2, z2, val1, val2);

    return intersection;
}

CIsoSurface.prototype.Interpolate = function(fX1, fY1, fZ1, fX2, fY2, fZ2, tVal1, tVal2)
{
    var interpolation = [-1,-1,-1];
    var mu;

    mu = (this.m_tIsoLevel - tVal1)/(tVal2 - tVal1);
    interpolation[0] = fX1 + mu*(fX2 - fX1);
    interpolation[1] = fY1 + mu*(fY2 - fY1);
    interpolation[2] = fZ1 + mu*(fZ2 - fZ1);

    return interpolation;
}

CIsoSurface.prototype.RenameVerticesAndTriangles = function()
{
    var nextID = 0;
    var i;
    var start = new Date().getTime();

    // Rename vertices.
    for (mapIterator in this.m_i2pt3idVertices) {
        this.m_i2pt3idVertices[mapIterator][3] = nextID;
        nextID++;
    }

    var end = new Date().getTime();
    var time = end - start;
    console.log("Time to RenameVerticesAndTriangles(1) : "+time*0.001);

    // Now rename triangles.
    for (vecIterator in this.m_trivecTriangles) {
        for (i = 0; i < 3; i++) {
            var newID = this.m_i2pt3idVertices[this.m_trivecTriangles[vecIterator][i]][3];
            this.m_trivecTriangles[vecIterator][i] = newID;
        }
    }

    end = new Date().getTime();
    var time = end - start;
    console.log("Time to RenameVerticesAndTriangles(2) : "+time*0.001);

    // Copy all the vertices and triangles into two arrays so that they
    // can be efficiently accessed.
    // Copy vertices.

    this.m_nVertices = nextID;
    this.m_ppt3dVertices = [];

    var i = 0;
    for (var mapIterator in this.m_i2pt3idVertices) {
        //var thisPoint = [-1,-1,-1];
        //this.m_ppt3dVertices.push(thisPoint);
        this.m_ppt3dVertices[i] = this.m_i2pt3idVertices[mapIterator];
        //this.m_ppt3dVertices[i][0] = this.m_i2pt3idVertices[mapIterator][0];
        //this.m_ppt3dVertices[i][1] = this.m_i2pt3idVertices[mapIterator][1];
        //this.m_ppt3dVertices[i][2] = this.m_i2pt3idVertices[mapIterator][2];
        i++;
    }

    end = new Date().getTime();
    var time = end - start;
    console.log("Time to RenameVerticesAndTriangles(3) : "+time*0.001);

    this.m_nTriangles = 0;
    this.m_piTriangleIndices = [];
    for (vecIterator in this.m_trivecTriangles) {
        this.m_piTriangleIndices.push(this.m_trivecTriangles[vecIterator][0]);
        this.m_piTriangleIndices.push(this.m_trivecTriangles[vecIterator][1]);
        this.m_piTriangleIndices.push(this.m_trivecTriangles[vecIterator][2]);
        this.m_nTriangles++;
    }

    end = new Date().getTime();
    var time = end - start;
    console.log("Time to RenameVerticesAndTriangles(4) : "+time*0.001);

    this.m_i2pt3idVertices = [];
    this.m_trivecTriangles = [];
    this.CalculateNormals();

    end = new Date().getTime();
    var time = end - start;
    console.log("Time to RenameVerticesAndTriangles(5) : "+time*0.001);

}

// debugging function
CIsoSurface.prototype.check_max_min_vertex_index_from_triangles = function() {

    var v_max = 0;
    var v1; 

    console.log("checking m_nTriangles=" + this.m_nTriangles + " triangles"); 
    console.log("         m_nVertices =" + this.m_nVertices + " vertices"); 

    for (var i = 0; i < this.m_nTriangles; i++) {

        // recall that m_piTriangleIndices is a list of vars
        // (indexing the vertices).

        v1 = this.m_piTriangleIndices[i];

        if (v1 > v_max) {
            v_max = v1;
        }

    }

    console.log("max vertex from triangle usage is: " + v_max);

}


// debugging function
CIsoSurface.prototype.check_max_min_vertices = function() {

    var max_x = -999, max_y = -999, max_z = -999; 
    var min_x =  999, min_y =  999, min_z =  999;

    console.log("checking m_nVertices=" + this.m_nVertices + " vertices"); 
    for (var i = 0; i < this.m_nVertices; i++) {

        if (this.m_ppt3dVertices[i][0] > max_x)
            max_x = this.m_ppt3dVertices[i][0]; 
        if (this.m_ppt3dVertices[i][1] > max_y)
            max_y = this.m_ppt3dVertices[i][1]; 
        if (this.m_ppt3dVertices[i][2] > max_z)
            max_z = this.m_ppt3dVertices[i][2]; 

        if (this.m_ppt3dVertices[i][0] < min_x)
            min_x = this.m_ppt3dVertices[i][0]; 
        if (this.m_ppt3dVertices[i][1] < min_y)
            min_y = this.m_ppt3dVertices[i][1]; 
        if (this.m_ppt3dVertices[i][2] < min_z)
            min_z = this.m_ppt3dVertices[i][2]; 

    }

    console.log("Debug: check_max_min_vertices (min and max x, y and z):" + min_x + " " << max_x);
    console.log(min_y + " " + max_y); 
    console.log(min_z + " " + max_z);

}


CIsoSurface.prototype.CalculateNormals = function()
{
    this.m_nNormals = this.m_nVertices;
    this.m_pvec3dNormals = [];

    //std::cout << "Calculating normals\n";

    var i;

    // Set all normals to 0.
    for (i = 0; i < this.m_nNormals; i++) {
        var normal = [-1,-1,-1];
        this.m_pvec3dNormals.push(normal);
        this.m_pvec3dNormals[i][0] = 0;
        this.m_pvec3dNormals[i][1] = 0;
        this.m_pvec3dNormals[i][2] = 0;
    }

    // Calculate normals.
    for (i = 0; i < this.m_nTriangles; i++) {
        var vec1 = [-1,-1,-1];
        var vec2 =[-1,-1,-1];
        var normal = [-1,-1,-1];
        var id0, id1, id2;
        id0 = this.m_piTriangleIndices[i*3];
        id1 = this.m_piTriangleIndices[i*3+1];
        id2 = this.m_piTriangleIndices[i*3+2];
        vec1[0] = this.m_ppt3dVertices[id1][0] - this.m_ppt3dVertices[id0][0];
        vec1[1] = this.m_ppt3dVertices[id1][1] - this.m_ppt3dVertices[id0][1];
        vec1[2] = this.m_ppt3dVertices[id1][2] - this.m_ppt3dVertices[id0][2];
        vec2[0] = this.m_ppt3dVertices[id2][0] - this.m_ppt3dVertices[id0][0];
        vec2[1] = this.m_ppt3dVertices[id2][1] - this.m_ppt3dVertices[id0][1];
        vec2[2] = this.m_ppt3dVertices[id2][2] - this.m_ppt3dVertices[id0][2];
        normal[0] = vec1[2]*vec2[1] - vec1[1]*vec2[2];
        normal[1] = vec1[0]*vec2[2] - vec1[2]*vec2[0];
        normal[2] = vec1[1]*vec2[0] - vec1[0]*vec2[1];
        this.m_pvec3dNormals[id0][0] += normal[0];
        this.m_pvec3dNormals[id0][1] += normal[1];
        this.m_pvec3dNormals[id0][2] += normal[2];
        this.m_pvec3dNormals[id1][0] += normal[0];
        this.m_pvec3dNormals[id1][1] += normal[1];
        this.m_pvec3dNormals[id1][2] += normal[2];
        this.m_pvec3dNormals[id2][0] += normal[0];
        this.m_pvec3dNormals[id2][1] += normal[1];
        this.m_pvec3dNormals[id2][2] += normal[2];
    }

    // Normalize normals.
    for (i = 0; i < this.m_nNormals; i++) {
        var length = Math.sqrt(this.m_pvec3dNormals[i][0]*this.m_pvec3dNormals[i][0] + this.m_pvec3dNormals[i][1]*this.m_pvec3dNormals[i][1] + this.m_pvec3dNormals[i][2]*this.m_pvec3dNormals[i][2]);
        this.m_pvec3dNormals[i][0] /= length;
        this.m_pvec3dNormals[i][1] /= length;
        this.m_pvec3dNormals[i][2] /= length;
    }
}

// PE
CIsoSurface.prototype.nTriangles = function() {

    return this.m_nTriangles;
}

/*
template <class T> var* CIsoSurface<T>::returnClippedIndices(const std::vector<Cartesian> &clip_positions_in, double clip_radius,  double x, double y, double z){

  if(clip_positions_in.size()==0) return returnIndices();
  if(m_nVertices<6) return returnIndices();

  std::vector<Cartesian> clip_positions = clip_positions_in;
  std::sort(clip_positions.begin(),clip_positions.end(),Cart_sort_x());

  double* vertices = returnVertices(x, y, z);

  var* indices = new var[m_nTriangles*3];
  var j=0;

  for (var i=0; i < m_nTriangles*3; i+=3) {
    Cartesian cart1(vertices+m_piTriangleIndices[i]);
    Cartesian cart2(vertices+m_piTriangleIndices[i+1]);
    Cartesian cart3(vertices+m_piTriangleIndices[i+2]);
    if(test_clip_pos(cart1,clip_positions,clip_radius)&&test_clip_pos(cart2,clip_positions,clip_radius)&&test_clip_pos(cart3,clip_positions,clip_radius)){
      indices[j++] = m_piTriangleIndices[i]; 
      indices[j++] = m_piTriangleIndices[i+2]; 
      indices[j++] = m_piTriangleIndices[i+1];
    }
  }

  m_nTriangles_clipped = j;

  return indices;
}
*/


CIsoSurface.prototype.returnIndices = function(){
    var indices = [];
    var j=0;

    for (var i=0; i < this.m_nTriangles*3; i+=3) {
        indices.push(this.m_piTriangleIndices[i]) 
            if(this.m_tIsoLevel>0){
                indices.push(this.m_piTriangleIndices[i+2]); 
                indices.push(this.m_piTriangleIndices[i+1]);
            } else {
                indices.push(this.m_piTriangleIndices[i+1]); 
                indices.push(this.m_piTriangleIndices[i+2]);
            }
    }

    return indices;
}

CIsoSurface.prototype.returnChickenIndices = function(){
    var indices = [];
    var idx = 0;
    for (var i=0; i < this.chicken_vertices.length; i+=3) {
        indices.push(0+idx);
        idx++;
    }
    return indices;
}

CIsoSurface.prototype.returnChickenNormals = function(){
    var norms = [];
    for (var i=0; i < this.chicken_vertices.length; i+=3) {
        norms.push(0.0);
        norms.push(0.0);
        norms.push(1.0);
    }
    return norms;
}

CIsoSurface.prototype.returnChickenVertices = function(){
    return this.chicken_vertices;
}

CIsoSurface.prototype.returnVertices = function(x, y, z, doChickenWire){

    //console.log(this.m_nVertices);
    //console.log(this.m_ppt3dVertices);

    var vertices = [];

    if(typeof(doChickenWire)!=="undefined" && doChickenWire){
        var chickLength = this.chicken_vertices.length;
        for (var i=0; i < chickLength; i+=3) {
            vertices.push(this.chicken_vertices[i]+x);
            vertices.push(this.chicken_vertices[i+1]+y);
            vertices.push(this.chicken_vertices[i+2]+z);
        }
        return vertices;
    }

    var j=0;
    if(this.cell){
        var vo = vec3.create();
        var nCellsX = this.cell.nx;
        var nCellsY = this.cell.ny;
        var nCellsZ = this.cell.nz;
        //console.log(this.cell);
        var RO = this.cell.matrix_orth;
        for (var i=0; i < this.m_nVertices; i++) {
            var v = vec3.create([(this.m_ppt3dVertices[i][0])/nCellsX,(this.m_ppt3dVertices[i][1])/nCellsY,(this.m_ppt3dVertices[i][2])/nCellsZ]);
            mat4.multiplyVec3(RO,v);
            vertices.push(v[0]+x);
            vertices.push(v[1]+y);
            vertices.push(v[2]+z);
            /*
               if(i==0){
               printMat(RO);
               console.log(v[0]+" "+v[1]+" "+v[2]);
               console.log("  "+(this.m_ppt3dVertices[i][0]+x)+" "+(this.m_ppt3dVertices[i][1]+y)+" "+(this.m_ppt3dVertices[i][2]+z));
               }
             */
        }
    } else {
        for (var i=0; i < this.m_nVertices; i++) {
            vertices.push(this.m_ppt3dVertices[i][0]+x);
            vertices.push(this.m_ppt3dVertices[i][1]+y);
            vertices.push(this.m_ppt3dVertices[i][2]+z);
        }

    }
    return vertices;
}

CIsoSurface.prototype.returnNormals_new = function(){

    var vertices = [];
    var j=0;
    for (var i=0; i < this.m_nVertices; i++) {
        vertices.push(Math.abs(this.m_tIsoLevel)/this.m_tIsoLevel*this.m_pvec3dNormals[i][0]);
        vertices.push(Math.abs(this.m_tIsoLevel)/this.m_tIsoLevel*this.m_pvec3dNormals[i][1]);
        vertices.push(Math.abs(this.m_tIsoLevel)/this.m_tIsoLevel*this.m_pvec3dNormals[i][2]);
    }
    return vertices;
}

CIsoSurface.prototype.returnTriangles = function(){

    var result = [];

    for (var i=0; i < this.m_nTriangles*3; i+=3) {

        var j   = this.m_piTriangleIndices[i]; 
        var jp  = this.m_piTriangleIndices[i+1]; 
        var jp2 = this.m_piTriangleIndices[i+2];

        var co1_c = this.m_ppt3dVertices[j];
        var co2_c = this.m_ppt3dVertices[jp];
        var co3_c = this.m_ppt3dVertices[jp2];

        var p1 = [];
        p1.push(co1_c);
        p1.push(co2_c);
        var p2 = [];
        p2.push(co1_c);
        p2.push(co3_c);
        var p3 = [];
        p3.push(co3_c);
        p3.push(co2_c);
        result.push(p1);
        result.push(p2);
        result.push(p3);

    }

    return result;
}

CIsoSurface.prototype.returnNormals = function(){
    /* 
     * This is an ugly hack to just make solid surface stuff
     * same as Paul's std::pair<Cartesian,Cartesian>  lines. The correct way is to
     * never have the std::pair<Cartesian,Cartesian>  vector at all, just use the
     * triangle data. However, it works.
     */

    var result = [];

    for (var i=0; i < this.m_nTriangles*3; i+=3) {

        var j   = this.m_piTriangleIndices[i]; 
        var jp  = this.m_piTriangleIndices[i+1]; 
        var jp2 = this.m_piTriangleIndices[i+2];

        var co1_c = this.m_pvec3dNormals[j];
        var co2_c = this.m_pvec3dNormals[jp];
        var co3_c = this.m_pvec3dNormals[jp2];

        var p1 = [];
        p1.push(co1_c);
        p1.push(co2_c);
        var p2 = [];
        p2.push(co1_c);
        p2.push(co3_c);
        var p3 = [];
        p3.push(co3_c);
        p3.push(co2_c);
        result.push(p1);
        result.push(p2);
        result.push(p3);

    }

    return result;
}
