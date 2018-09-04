#include "platform_unix.inl"

/*
struct DirectoryWatch {
    FSEventStreamRef fs_stream;
    Url path;
    DynamicArray<Url> files;
};

static DynamicArray<DirectoryWatch> s_watched_dirs;

bool AddDirectoryWatch(CString path_to_watch, void (*callback)(FileEvent)) {
    //Define variables and create a CFArray object containing
     CFString objects containing paths to watch.
     




    CFStringRef path = CFStringCreateWithCString(NULL, path_to_watch, CFStringBuiltInEncodings::kCFStringEncodingUTF8);
    CFArrayRef path_array = CFArrayCreate(NULL, (const void **)&path, 1, NULL);
    FSEventStreamRef stream;
    const CFAbsoluteTime latency = 1.0; // Latency in seconds
    




    // Create the stream, passing in a callback
    stream = FSEventStreamCreate(NULL,
                                 &mycallback,
                                 NULL,
                                 path_array,
                                 kFSEventStreamEventIdSinceNow,
                                 latency,
                                 kFSEventStreamCreateFlagNone);
    




    FSEventStreamScheduleWithRunLoop(stream, CFRunLoopGetCurrent(), kCFRunLoopDefaultMode);
    FSEventStreamStart(stream);
    




    return true;
}

bool RemoveWatch(CString path_to_watch) {
    CString dir = get_directory(path_to_watch);
    CString file = get_file(path_to_watch);
    




    FSEventStreamStop(stream);
    FSEventStreamInvalidate(stream);
    FSEventStreamRelease(stream);
}
*/
