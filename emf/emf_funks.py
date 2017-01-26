from . import os, np, pd

import emf_class

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
    if('.' in os.path.basename(filename_if_needed)):
        filename_if_needed = filename_if_needed[:filename_if_needed.rfind('.')]
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
            if('.' in os.path.basename(path)):
                path = path[:path.rfind('.')]
            #put intended extension on the end and return
            return(path + extension)
    else: #otherwise just use filename_if_needed and extension
        return(filename_if_needed + extension)

def _check_extension(file_path, correct_ext, message):
    """Check that a file path has the desired extension, raising an error if
    not and appending the correct extension if no extension is present.
    args:
        file_path - string, a target file path
        correct_ext - string, the correct extension for the target path,
                      with or without the period
        message - string, error message if the extention is wrong
    returns:
        file_path - string, valid file path"""

    #manage period in correct_ext
    if(not (correct_ext[0] == '.')):
        correct_ext = '.' + correct_ext

    #if it has an extension, check it
    if('.' in os.path.basename(file_path)):
        #check that the file exists
        if(not os.path.isfile(file_path)):
            raise(emf_class.EMFError("""
            The file path:
                "%s"
            is not recognized as an existing file.""" % file_path))
        #get extension and check it against correct_ext
        ext = file_path[file_path.rfind('.'):]
        if(ext != correct_ext):
            raise(emf_class.EMFError(message))
        else:
            return(file_path)
    #no extension, just append the correct extension and check file exists
    else:
        file_path += correct_ext
        #check that the file exists
        if(not os.path.isfile(file_path)):
            raise(emf_class.EMFError("""
            The file path:
                "%s"
            is not recognized as an existing file.""" % file_path))
        return(file_path)

def _is_number(s):
    """Check if an element can be converted to a float, returning True
    if it can and False if it can't"""
    if((s is False) or (s is True)):
        return(False)
    try:
        float(s)
    except(ValueError, TypeError):
        return(False)
    else:
        return(True)

def _is_int(x):
    if(_is_number(x)):
        if(float(x) == int(float(x))):
            return(True)
    return(False)

def _sig_figs(v, figs):
    if(v == 0):
        return(0)
    if(pd.isnull(v)):
        return(v)
    w = round(v, int(figs - np.ceil(np.log10(np.abs(v)))))
    if(w == round(w)):
        return(int(w))
    return(w)

def _check_intable(f):
    """If a float is a whole number, convert it to an integer"""
    if(_is_int(f)):
        return(int(f))
    else:
        return(float(f))

#Flatten a list of lists. Only works on a list and the first level of
#sublists, not recursively on all levels
_flatten = lambda L: [item for sublist in L for item in sublist]

_to_alphanumeric = lambda s: ''.join(ch for ch in s if ch.isalnum())

def _path_str_condition(s):
	s = s.replace(' ', '-')
	return(_to_alphanumeric(s))

def _Levenshtein_group(V, W):
    """Match each string in vector V to a string in vector W, using global
    least Levenshtein distance for each match.
        V - vector of strings
        W - vector of strings, can be longer than V but not shorter
    returns:
        list of matches in the same order as V, comprising elements of W"""

    #set up some loop variables
    vmask = np.ones((len(V),), dtype = bool)
    vrange = np.arange(len(V))
    wmask = np.ones((len(W),), dtype = bool)
    wrange = np.arange(len(W))
    matches = ['']*len(V)
    #calculate distance between all pairs of words
    M = np.zeros((len(V),len(W)))
    for i in vrange:
        for j in wrange:
            M[i,j] = _Levenshtein_distance(V[i],W[j])
    #find a match in W for each word in V with minimum global distances
    for n in vrange:
        vr = vrange[vmask]
        wr = wrange[wmask]
        #find the minimum distance of any pairs left
        i_min = vr[0]
        j_min = wr[0]
        min_dist = M[i_min,j_min]
        for i in vr:
            for j in wr:
                if(M[i,j] < min_dist):
                    min_dist = M[i,j]
                    i_min = i
                    j_min = j
        #store the match and update the masks
        matches[i_min] = W[j_min]
        vmask[i_min] = False
        wmask[j_min] = False

    return(matches)

def _Levenshtein_find(s, V):
    """find the string in V that most clostly resembles s, as measured by
    the Levenshtein distance
    args:
        s - string
        V - iterable of strings
    returns:
        string in V most similar to s
        index of that string in V"""
    d = np.zeros((len(V),))
    for (i,t) in enumerate(V):
        d[i] = _Levenshtein_distance(s,t)
    idx = np.argmin(d)
    return(V[idx], idx)

def _Levenshtein_distance(s, t):
    """Calculate string similarity metric with a basic edit distance
    routine to find the Levenshtein Distance, without case sensitivity
    args:
        s - string 1
        t - string 2
    returns:
        d - distance metric, count of edit operations required
    information source:
        http://web.archive.org/web/20081224234350/http://www.dcs.shef.ac.uk/~sam/stringmetrics.html#Levenshtein"""
    #check empty strings
    if(len(s) == 0):
        return(len(t))
    elif(len(t) == 0):
        return(len(s))
    #lower case
    s = s.lower()
    t = t.lower()
    #initialize grid
    ls = len(s)
    lt = len(t)
    D = np.zeros((ls,lt))
    if(s[0] != t[0]):
        D[0,0] = 1.0
    D[:,0] = np.arange(D[0,0], ls + D[0,0])
    D[0,:] = np.arange(D[0,0], lt + D[0,0])
    #vector to store edit operation scores
    e = np.zeros((3,))
    for i in range(1,ls):
        for j in range(1,lt):
            e[0] = D[i-1,j-1]
            if(s[i] != t[j]):
                e[0] += 1
            e[1] = D[i-1,j] + 1
            e[2] = D[i,j-1] + 1
            D[i,j] = np.min(e)
    return(D[-1,-1])
