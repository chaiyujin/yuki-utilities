#pragma once

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

#include <thread>
#include <mutex>
#include <condition_variable>
#include <queue>

namespace yuki
{
	#ifndef __YUKI_UTILS_LOG_H__
	#define __YUKI_UTILS_LOG_H__

	// #define CHECK_LEVEL_NONE // uncomment this line to disable runtime check


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
}

namespace yuki
{
	template <class T>
	class SafeQueue
	{
	public:
		SafeQueue(void): q_(), m_(), c_(), size_(0) {}
		~SafeQueue(void){}
		// Add an element to the queue.
		void push(T t)
		{
			std::lock_guard<std::mutex> lock(m_);
			q_.push(t);
			size_++;
			c_.notify_one();
		}
		// Get the "front"-element.
		// If the queue is empty, wait till a element is avaiable.
		void pop(void)
		{
			std::unique_lock<std::mutex> lock(m_);
			while (q_.empty())
			{
				// release lock as long as the wait and reaquire it afterwards.
				c_.wait(lock);
			}
			q_.pop();
			size_--;
		}
		T &front()
		{
			std::unique_lock<std::mutex> lock(m_);
			while (q_.empty())
			{
				// release lock as long as the wait and reaquire it afterwards.
				c_.wait(lock);
			}
			return q_.front();
		}

		size_t size() const { return size_; }

	private:
		std::queue<T> q_;
		mutable std::mutex m_;
		std::condition_variable c_;
		size_t size_;
	};
}


#define NAMESPACE_BEGIN(name) namespace name {
#define NAMESPACE_END(name) }