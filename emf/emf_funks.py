import os

def _path_manage(filename_if_needed, extension, **kwargs):
    """This function takes a path string through the kwarg 'path' and
    returns a path string with a file name at it's end, to save a file
    at that location. If the path string is a directory (ends with a slash
    or doesn't but has preceding directory elements) a new path string is
    returned with the 'filename_if_needed' and 'extension' arguments
    appended. If the path string already has a file name at its end, the
    input extension will replace any preexisting one. If no path string is
    passed in via the keyword argument 'path', the returned path is simply
    the input filename_if_needed with the input extension at its end. If the
    last part of the path string is supposed to be the filename, include an
    extension.
    args:
        filename_if_needed - name of file if not in path string
        extension - file extension
    kwargs:
        path - string, destination/filename for saved file(s)
    returns:
        conditioned path string"""
    #remove extensions from filename_if_needed
    if('.' in filename_if_needed):
        filename_if_needed = filename_if_needed[:filename_if_needed.index('.')]
    #make sure the extension has a period at its beginning
    if(extension):
        if(extension[0] != '.'):
            extension = '.' + extension
    #check if there is a path kwarg
    if('path' in kwargs):
        path = kwargs['path']
        #check if the input path is a directory
        if(os.path.isdir(path)):
            #join the directory string, file string, and extension string
            return(os.path.join(path, filename_if_needed + extension))
        else:
            #path string interpreted as path to a file, not directory
            #strip any extensions
            if('.' in path):
                path = path[:path.rfind('.')]
            #put intended extension on the end and return
            return(path + extension)
    else: #otherwise just use filename_if_needed and extension
        return(filename_if_needed + extension)

def _check_extension(file_path, correct_ext, message):
    """Check that a file path has the desired extension, raising an error if
    not and appending the correct extension if no extension is present.
    args:
        file_path - a target file path
        correct_ext - the correct extension for the target path, with or
                      without the period
        message - error message if the extention is wrong
    returns:
        file_path"""
    if(correct_ext[0] == '.'):
        correct_ext = correct_ext[1:]
    if('.' in file_path):
        ext = file_path[file_path.index('.')+1:]
        if(not (correct_ext == ext)):
            raise(emf_class.EMFError(message))
    else:
        file_path += '.' + correct_ext
    return(file_path)

def _is_number(s):
    """Check if an element can be converted to a float, returning `True`
    if it can and `False` if it can't"""
    if((s is False) or (s is True)):
        return(False)
    try:
        float(s)
    except(ValueError, TypeError):
        return(False)
    else:
        return(True)

def _check_intable(f):
    """If a float is a whole number, convert it to an integer"""
    if(float(f) == int(f)):
        return(int(f))
    else:
        return(float(f))

#Flatten a list of lists. Only works on a list and the first level of
#sublists, not recursively on all levels
_flatten = lambda L: [item for sublist in L for item in sublist]
