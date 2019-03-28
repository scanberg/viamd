/*
 * basE91 command line front-end
 *
 * Copyright (c) 2000-2006 Joachim Henke
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *  - Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *  - Neither the name of Joachim Henke nor the names of his contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <fcntl.h>
#include <unistd.h>
#endif
#include <getopt.h>
#include "base91.h"

#define FLG_D 1
#define FLG_V 2
#define FLG_VV 4

static char status[32];
static const char *progname;
static char *ibuf, *obuf;
static size_t ibuf_size, llen;
static struct basE91 b91;

static void stream_b91enc_p(void)
{
	size_t itotal = 0;
	size_t ototal = 0;
	size_t s;

	while ((s = fread(ibuf, 1, ibuf_size, stdin)) > 0) {
		itotal += s;
		s = basE91_encode(&b91, ibuf, s, obuf);
		ototal += s;
		fwrite(obuf, 1, s, stdout);
	}
	s = basE91_encode_end(&b91, obuf);	/* empty bit queue */
	ototal += s;
	fwrite(obuf, 1, s, stdout);

	sprintf(status, "\t%.2f%%\n", itotal ? (float) ototal / itotal * 100.0 : 1.0);
}

static void stream_b91enc_w(void)
{
	size_t l = llen;
	size_t ltotal = 0;
	size_t i, s;
	char x;

	while ((s = fread(ibuf, 1, ibuf_size, stdin)) > 0) {
		s = basE91_encode(&b91, ibuf, s, obuf);
		for (i = 0; l <= s; l += llen) {
			x = obuf[l];
			obuf[l] = '\0';
			puts(obuf + i);
			++ltotal;
			obuf[l] = x;
			i = l;
		}
		fwrite(obuf + i, 1, s - i, stdout);
		l -= s;
	}
	s = basE91_encode_end(&b91, obuf);
	if (s || l < llen) {
		obuf[s] = '\0';
		if (s > l) {
			x = obuf[1];
			obuf[1] = '\0';
			puts(obuf);
			++ltotal;
			obuf[0] = x;
		}
		puts(obuf);
		++ltotal;
	}

	sprintf(status, "\t%lu lines\n", (unsigned long) ltotal);
}

static void stream_b91dec(void)
{
	size_t s;

	while ((s = fread(ibuf, 1, ibuf_size, stdin)) > 0) {
		s = basE91_decode(&b91, ibuf, s, obuf);
		fwrite(obuf, 1, s, stdout);
	}
	s = basE91_decode_end(&b91, obuf);	/* empty bit queue */
	fwrite(obuf, 1, s, stdout);

	sprintf(status, "done\n");
}

static int init_flags(const char *p)
{
	size_t l = strlen(p);

	if (l > 5) {
		progname = p + l - 6;
		if (!strcmp(progname, "b91enc"))
			return 0;
		if (!strcmp(progname, "b91dec"))
			return FLG_D;
	}
	llen = 76;
	progname = "base91";

	return 0;
}

int main(int argc, char **argv)
{
	size_t buf_size = 65536;	/* buffer memory defaults to 64 KiB */
	int flags = init_flags(*argv);
	char *ifile = "from standard input";
	char *ofile = NULL;
	int opt;
	struct option longopts[8] = {
		{"decode", no_argument, NULL, 'd'},
		{"output", required_argument, NULL, 'o'},
		{"verbose", no_argument, NULL, 'v'},
		{"wrap", required_argument, NULL, 'w'},
		{"help", no_argument, NULL, 'h'},
		{"version", no_argument, NULL, 'V'},
		{NULL, 0, NULL, 0}
	};

	while ((opt = getopt_long(argc, argv, "dem:o:vw:hV", longopts, NULL)) != -1)
		switch (opt) {
		case 'd':
			flags |= FLG_D;
			break;
		case 'e':
			flags &= ~FLG_D;
			break;
		case 'm':
			{
				char *t;
				long l = strtol(optarg, &t, 0);

				if (t == optarg || strlen(t) > 1 || l < 0) {
					fprintf(stderr, "invalid SIZE argument: `%s'\n", optarg);
					return EXIT_FAILURE;
				}
				buf_size = l;
				switch (*t | 32) {
				case ' ':
				case 'b':
					break;
				case 'k':
					buf_size <<= 10;
					break;
				case 'm':
					buf_size <<= 20;
					break;
				default:
					fprintf(stderr, "invalid SIZE suffix: `%s'\n", t);
					return EXIT_FAILURE;
				}
			}
			break;
		case 'o':
			if (strcmp(optarg, "-"))
				ofile = optarg;
			break;
		case 'v':
			flags |= (flags & FLG_V) ? FLG_VV : FLG_V;
			break;
		case 'w':
			{
				char *t;
				long l = strtol(optarg, &t, 0);

				if (*t || l < 0) {
					fprintf(stderr, "invalid number of columns: `%s'\n", optarg);
					return EXIT_FAILURE;
				}
				llen = l;
			}
			break;
		case 'h':
			printf("Usage: %s [OPTION]... [FILE]\n"
				"basE91 encode or decode FILE, or standard input, to standard output.\n", progname);
			puts("\n  -d, --decode\t\tdecode data\n"
				"  -m SIZE\t\tuse SIZE bytes of memory for buffers (suffixes b, K, M)\n"
				"  -o, --output=FILE\twrite to FILE instead of standard output\n"
				"  -v, --verbose\t\tverbose mode\n"
				"  -w, --wrap=COLS\twrap encoded lines after COLS characters (default 76)\n"
				"  --help\t\tdisplay this help and exit\n"
				"  --version\t\toutput version information and exit\n\n"
				"With no FILE, or when FILE is -, read standard input.");
			return EXIT_SUCCESS;
		case 'V':
			printf("%s 0.6.0\nCopyright (c) 2000-2006 Joachim Henke\n", progname);
			return EXIT_SUCCESS;
		default:
			fprintf(stderr, "Try `%s --help' for more information.\n", *argv);
			return EXIT_FAILURE;
		}

	if (flags & FLG_D) {
		ibuf_size = (buf_size - 1) << 3;
		if (ibuf_size < 15) {
			fputs("SIZE must be >= 3 for decoding\n", stderr);
			return EXIT_FAILURE;
		}
		ibuf_size /= 15;
	} else {
		ibuf_size = (buf_size - 2) << 4;
		if (ibuf_size < 29) {
			fputs("SIZE must be >= 4 for encoding\n", stderr);
			return EXIT_FAILURE;
		}
		ibuf_size /= 29;
	}

	if (optind < argc && strcmp(argv[optind], "-")) {
		ifile = argv[optind];
		if (freopen(ifile, "r", stdin) != stdin) {
			perror(ifile);
			return EXIT_FAILURE;
		}
	}
	if (ofile)
		if (freopen(ofile, "w", stdout) != stdout) {
			perror(ofile);
			return EXIT_FAILURE;
		}

	if (flags & FLG_VV)
		fprintf(stderr, "using %lu bytes for buffers; input buffer: %lu bytes\n", (unsigned long) buf_size, (unsigned long) ibuf_size);
	obuf = malloc(buf_size);
	if (!obuf) {
		fputs("failed to allocate buffer memory\n", stderr);
		return EXIT_FAILURE;
	}

	basE91_init(&b91);
#ifdef _WIN32
	_setmode(_fileno(stdin), _O_BINARY);
#endif

	if (flags & FLG_D) {
#ifdef _WIN32
		_setmode(_fileno(stdout), _O_BINARY);
#endif
		ibuf = obuf + 1;	/* create overlapping buffers to use memory efficiently */
		if (flags & FLG_V)
			fprintf(stderr, "decoding %s ...", ifile);
		stream_b91dec();
	} else {
		ibuf = obuf + buf_size - ibuf_size;	/* partial overlap */
		if (flags & FLG_V)
			fprintf(stderr, "encoding %s ...", ifile);
		if (llen)
			stream_b91enc_w();
		else
			stream_b91enc_p();
	}
	free(obuf);

	if (flags & FLG_V)
		fputs(status, stderr);

	return EXIT_SUCCESS;
}
