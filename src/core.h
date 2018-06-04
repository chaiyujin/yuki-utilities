#pragma once

#define NAMESPACE_BEGIN(name) namespace name {
#define NAMESPACE_END(name) }

#ifndef __YUKI_UTILS_LOG_H__
#define __YUKI_UTILS_LOG_H__

// #define CHECK_LEVEL_NONE // uncomment this line to disable runtime check

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#define YUKI_LOG(fmt, ...) { printf("[LOG] "); printf(fmt, ##__VA_ARGS__); }

#ifdef NDEBUG
#define YUKI_DEBUG(fmt, ...) {}
#else
#define YUKI_DEBUG(fmt, ...) { printf("[DEBUG] "); printf(fmt, ##__VA_ARGS__); }
#endif

#define YUKI_ERROR(fmt, ...) {\
	printf("Error in [%s] (line: %d) : ", __FUNCTION__, __LINE__);\
	printf(fmt, ##__VA_ARGS__); }
#define YUKI_ERROR_EXIT(fmt, ...) {\
	printf("Error in [%s] (line: %d) : ", __FUNCTION__, __LINE__);\
	printf(fmt, ##__VA_ARGS__);\
	getchar(); exit(1); }
#define CHECK(assertion) {\
		bool flag = assertion;\
		if (!flag) {\
			printf("Error in [%s] (line: %d) : CHECK fail.\n", __FUNCTION__, __LINE__);\
			getchar(); exit(1);\
		}\
	}

#endif // !__YUKI_UTILS_LOG_H__