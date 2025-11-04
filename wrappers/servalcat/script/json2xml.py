def json2xml(json_obj, tag_name=None, tag_name_subroot="subroot"):
# https://gist.github.com/bewestphal/0d2f884b3327f31b86122b5f6a38b3e2
    def tag_name_clear(tag_name):
        tag_name_clean = tag_name.replace(" ", "_")
        tag_name_clean = tag_name_clean.replace("-", "minus")
        tag_name_clean = tag_name_clean.replace(",", "")
        tag_name_clean = tag_name_clean.replace(".", "")
        tag_name_clean = tag_name_clean.replace("(", "")
        tag_name_clean = tag_name_clean.replace(")", "")
        tag_name_clean = tag_name_clean.replace("/", "")
        tag_name_clean = tag_name_clean.replace("|", "")
        tag_name_clean = tag_name_clean.replace("*", "")
        tag_name_clean = tag_name_clean.replace(">", "_gt_")
        tag_name_clean = tag_name_clean.replace("<", "_lt_")
        tag_name_clean = tag_name_clean.replace("=", "_eq_")
        tag_name_clean = tag_name_clean.replace(":", "_")
        tag_name_clean = tag_name_clean.replace(";", "_")
        return tag_name_clean
    result_list = list()
    json_obj_type = type(json_obj)
    if tag_name:
        tag_name_clean = tag_name_clear(tag_name)
    else:
        tag_name_clean = tag_name_subroot
    if json_obj_type is list:
        for sub_elem in json_obj:
            result_list.append("\n<%s>" % (tag_name_clean))
            result_list.append(json2xml(sub_elem, tag_name=tag_name))
            result_list.append("</%s>" % (tag_name_clean))
        return "".join(result_list)
    if json_obj_type is dict:
        for tag_name in json_obj:
            if tag_name: tag_name_clean = tag_name_clear(tag_name)
            sub_obj = json_obj[tag_name]
            if isinstance(sub_obj, list):
                result_list.append(json2xml(sub_obj, tag_name=tag_name))
            elif isinstance(sub_obj, dict):
                result_list.append("\n<%s>" % (tag_name_clean))
                result_list.append(json2xml(sub_obj, tag_name=tag_name))
                result_list.append("\n</%s>" % (tag_name_clean))
            else:
                result_list.append("\n<%s>" % (tag_name_clean))
                result_list.append(json2xml(sub_obj, tag_name=tag_name))
                result_list.append("</%s>" % (tag_name_clean))
        return "".join(result_list)
    return "%s" % json_obj
