import json


def write(data, filename):
    '''
    metadata is plain text json file, which can be viewed using:
    http://jsonviewer.stack.hu/

    {
    "type":"block",
    "data":[
        {
            "file":"01/",
            "pos":[0, 0, 0]
        },
        {
            "file":"02/",
            "pos":[0, 10, 0]
        },
        {
            "file":"03/",
            "pos":[0, 30, 0]
        }
    ]
    }
    '''
    with open(filename, 'w') as outfile:
        json.dump(data, outfile, indent=4)


# get data block from file
def parse_data(filename, jsontype=None):
    content = parse_all(filename, jsontype)
    if content is not None:
        rtn = _finalize(content['data'], jsontype)
        return rtn
    else:
        return None


# get all
def parse_all(filename, jsontype=None):
    content = None
    with open(filename) as fp:
        content = json.load(fp)
    if content is not None:
        # handle exception
        if jsontype is not None and content['type'] != jsontype:
            print('[json helper]\tjson file type must be', jsontype)
            return None
        return content
    else:
        return None


def _finalize(rawjson, jsontype):
    # so far, data conversion is handled in their own position.
    if jsontype == "split":
        return rawjson
    elif jsontype == "DiMorSC":
        return rawjson
    elif jsontype == "block":
        return rawjson
    else:
        print("[json helper]\tUnsupported json type")
        return None


if __name__  == "__main__":
    pass
