from collections import Counter



# Your hits list with 25 IDs
id_list = [
    "FBgn0053207",
    "FBgn0010341",
    "FBgn0004580",
    "FBgn0032487",
    "FBgn0001250",
    "FBgn0011656",
    "FBgn0005771",
    "FBgn0038056",
    "FBgn0038056",
    "FBgn0038056",
    "FBgn0038056",
    "FBgn0038056",
    "FBgn0038056",
    "FBgn0030789",
    "FBgn0024321",
    "FBgn0024321",
    "FBgn0016797",
    "FBgn0037816",
    "FBgn0037816",
    "FBgn0004606",
    "FBgn0002931"
]

from collections import Counter



# Count the occurrences of each ID
id_counts = Counter(id_list)

# Initialize counters
appear_once = 0
appear_more_than_once = 0

# Iterate through the counts
for id, count in id_counts.items():
    if count == 1:
        appear_once += 1
    elif count > 1:
        appear_more_than_once += 1

print("IDs that appear only once:", appear_once)
print("IDs that appear more than once:", appear_more_than_once)


# Print the counts for each category
print("IDs that appear 0 times:", 592-(appear_once+appear_more_than_once), "{:.3f}".format((592-(appear_once+appear_more_than_once))/592))
print("IDs that appear 1 time:", appear_once, "{:.3f}".format(appear_once/592))
print("IDs that appear more than 2 times:", appear_more_than_once, "{:.3f}".format(appear_more_than_once/592))
