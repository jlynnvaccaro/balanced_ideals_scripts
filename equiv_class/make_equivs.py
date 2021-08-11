import equivalence_class
import utils
import time

# for i in range(1,9):


A = equivalence_class.equiv_tree(1)
start = time.time()
for i in range(1,8):
    A.up_cartan(i)
    A.make_equiv_class()
    end = time.time()
    print("Time diff {}: {:.4f}s".format(i,end-start))
    start = end
    # d = {}
    # utils.add_header(d)
    # cartan = "A"+str(i)
    # filename = cartan + ".json"
    # A = equivalence_class.equiv_tree(i)
    # A.make_equiv_class()
    # d["cartan_type"] = cartan
    # d["equiv_class"] = equivalence_class.to_str_list(A.equiv_nodes())
    # d["num_equiv_class"] = len(d["equiv_class"])
    # utils.save_json_data(filename,d)

start = time.time()
for i in range(1,8):
    A = equivalence_class.equiv_tree(i)
    A.make_equiv_class()
    end = time.time()
    print("{}: {:.4f}s".format(i,end-start))
    start = end