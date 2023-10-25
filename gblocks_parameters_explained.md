# Gblocks - parameters explained

For more details refer to: https://home.cc.umanitoba.ca/~psgendb/doc/Castresana/Gblocks_documentation.html
This is a simplify version related to our use.

## **1. `b1` (Minimum Number of Sequences for a Conserved Position):**

This parameter refers to the minimum number of sequences that need to have a conserved position to consider it as such in the alignment. Here, "conserved" implies that a specific position in the alignment is occupied by the same (or very similar) nucleotide or amino acid across multiple sequences.

For example, consider you have an alignment of 10 sequences. If you set **`-b1=7`**, it means that at least 7 out of the 10 sequences should have the same nucleotide or amino acid at a specific position for that position to be deemed conserved. This parameter ensures that only those positions in the alignment that have a significant level of conservation across the sequences are retained.

The default value for this parameter is usually half the total number of sequences plus one. For instance, with 10 sequences, the default would be 6 (10/2 + 1).

To understand the importance of this, consider the following simplified alignment example with 10 sequences:

```makefile

Seq1:  A-TA
Seq2:  ACTA
Seq3:  A-TA
Seq4:  AATA
Seq5:  AGTA
Seq6:  A-TA
Seq7:  A-TA
Seq8:  A-TA
Seq9:  ACTA
Seq10: A-TA

```

For position 2, we have the following nucleotides: **`-CCCCGT---C`**. If **`-b1`** is set to 7, then this position is considered conserved because 7 sequences have a **`C`** (even if 3 sequences have a gap or another nucleotide).

In choosing the value for this parameter, one must consider the trade-off between strictness (high conservation across most sequences) and inclusivity (allowing positions that might be conserved in a smaller subset of sequences).

### From practical standpoint:

The higher the number, the more strict the limit for a position to be designated as conserved.

To put it more clearly:

- A **higher value** for **`b1`** means that a larger number of sequences must have the same nucleotide or amino acid at a particular position for that position to be considered conserved. This makes the criteria for designating a position as "conserved" stricter.
- Conversely, a **lower value** for **`b1`** means fewer sequences need to have the same nucleotide or amino acid for a position to be considered conserved. This is a more lenient setting.

For instance, in a dataset of 100 sequences:

- If **`b1`** is set to 90, then at least 90 of those sequences need to have the same nucleotide or amino acid at a position for it to be considered conserved.
- If **`b1`** is set to 50, then only 50 of those sequences need to have the same nucleotide or amino acid at that position for it to be seen as conserved.

So, to reiterate, a higher **`-b1`** value makes the criteria stricter, while a lower **`-b1`** value makes it more lenient.

### Limits

Gblocks has some built-in limitations and default values to ensure the results remain biologically meaningful.

For the **`-b1`** parameter (Minimum Number Of Sequences For A Conserved Position):

- The **minimum value** allowed is half the total number of sequences (rounded down if the number of sequences is odd). This ensures that a position cannot be designated as "conserved" if it's conserved in less than half the sequences.
- The **maximum value** is the total number of sequences minus one. This means that, at its strictest setting, all sequences except one need to have the same nucleotide or amino acid at a position for it to be considered conserved.
- The **default value**, if you don't specify anything, is set to 85% of the total number of sequences (rounded down). This offers a balanced setting that's neither too lenient nor too strict.

To put it into perspective, if you have 100 sequences:

- The minimum value for **`b1`** is 50.
- The maximum value for **`b1`** is 99.
- The default value for **`b1`** is 85 (if not specified).

## 02. **`-b2`: Minimum Number Of Sequences For A Flanking Position**

This parameter determines how many sequences must have a conserved position in order for that position to serve as a "flank" to a block. Remember, blocks in Gblocks are the contiguous stretches of conserved positions that pass the criteria. These blocks are flanked by conserved positions which serve as their boundaries.

Essentially, **`-b2`** specifies the minimum number of sequences that must have a conserved position in order for that position to serve as the start or end of a block. The rationale is that these flanking positions give context and boundaries to the blocks and should thus be reasonably conserved.

For example, if you set **`-b2`** to 100, this means that for a position to serve as a flanking position, it must be conserved in at least 100 sequences of your alignment.

Again, like **`-b1`**, the higher you set **`-b2`**, the more stringent or strict the criteria for flanking positions. The opposite is also true: a lower **`-b2`** value makes the criteria for designating flanking positions more relaxed.

### Limits

The **`-b2`** parameter, which specifies the Minimum Number Of Sequences For A Flanking Position, has a lower limit. This lower limit is equal to the number you've set for the **`-b1`** parameter (Minimum Number Of Sequences For A Conserved Position). This is because a flanking position, by definition, should also be conserved, so it can't be less conserved than the general conserved positions.

In simpler terms, you can't have a flanking position that is less conserved than any random conserved position in your alignment.

Therefore:

- The lower limit for **`b2`** is the value you set for **`b1`**.
- The upper limit is the total number of sequences in your alignment.

If you're trying to set Gblocks parameters in a relaxed manner to retain more of your alignment, you would typically set **`-b2`** to a value close to the value of **`-b1`**, which would be around half the number of sequences in the alignment (as discussed earlier for a relaxed approach).

For the **`-b2`** parameter in Gblocks, which specifies the Minimum Number Of Sequences For A Flanking Position, the default value is set to 85% of the number of sequences in your alignment.

So, if you have 100 sequences in your alignment, the default setting for **`-b2`** would consider a position as a flanking position if it is present in at least 85 sequences.

## 03. **`-b3` (Maximum Number Of Contiguous Nonconserved Positions):**

This parameter determines how many non-conserved positions (in a row) are allowed in the final blocks. Gblocks considers positions in a sequence alignment to be "conserved," "flanking," or "non-conserved." After identifying conserved and flanking positions based on the previous parameters (**`-b1`** and **`-b2`**), Gblocks will then look for stretches of non-conserved positions. If the number of these contiguous non-conserved positions exceeds the limit set by **`-b3`**, they will be excluded from the final blocks.

In simpler terms, it allows you to specify the maximum length of a gap or "hole" in the alignment that you're willing to tolerate. This helps to eliminate regions of the alignment that are too variable to be reliably informative.

For example, if you set **`-b3=8`**, then Gblocks will allow up to 8 contiguous non-conserved positions in a block. If there's a stretch of 9 or more non-conserved positions in a row, that entire stretch will be removed from the final block.

**Default Value:** The default setting for **`-b3`** is 8, meaning Gblocks will tolerate up to 8 contiguous non-conserved positions by default.

It's important to note that setting this value too high can include regions of the alignment that are too variable, whereas setting it too low can exclude potentially informative parts of the alignment. As with other parameters, it's a balance between stringency and keeping potentially informative regions. Adjusting this parameter (along with the others) and checking the impact on your results can help you find an optimal setting for your specific dataset.

## 04 **`-b4` (Minimum Length Of A Block):**

This parameter sets the minimum number of contiguous positions that a block must have to be kept in the final alignment. In other words, once Gblocks identifies conserved and flanking regions and excludes stretches of non-conserved positions based on the **`-b3`** setting, it will further ensure that each resulting block has at least as many positions as specified by **`-b4`**.

The idea behind this is to prevent very short, potentially spurious blocks from being included in the final alignment. Short blocks might not be phylogenetically informative or could be the result of random sequence similarities, especially in large alignments.

For instance, if **`-b4=10`**, then any block that's shorter than 10 positions will be excluded from the final alignment.

**Default Value:** The default value for **`-b4`** is set to 10. So, by default, Gblocks requires blocks to be at least 10 positions long.

When adjusting this parameter, it's crucial to consider the nature of your dataset and the potential biological significance of shorter blocks. For some datasets, allowing shorter blocks might be beneficial if those regions are known to be of particular importance or interest. Conversely, in other cases, shorter blocks might just introduce noise into the analysis, so a stricter **`-b4`** setting might be preferred.

## 05 **`-b5` (Allowed Gap Positions):**

This parameter specifies how gaps in the alignment are treated when determining conserved positions. Remember, in the context of multiple sequence alignments, gaps are usually represented by the **`-`** character and signify that there is no corresponding residue (like a nucleotide or amino acid) in that position for a particular sequence.

Here's a breakdown of the possible settings for **`-b5`**:

1. **b5=a (All gaps allowed)**:
    - This is the most relaxed setting regarding gaps.
    - Gblocks will allow columns (positions) in the alignment that consist only of gaps.
    - This might be beneficial in cases where gaps themselves are informative or when using methods that can accommodate entirely gapped columns.
2. **b5=h (Half)**:
    - This is the default setting.
    - Positions with gaps in more than half the number of sequences are not allowed. In other words, if more than 50% of the sequences have a gap at a particular position, that position will not be considered conserved, and it might be excluded depending on the other parameters.
    - This setting is a balanced approach, ensuring that a majority of sequences have data for a given position while still allowing some flexibility.
3. **b5=n (No gaps allowed)**:
    - This is the most strict setting regarding gaps.
    - Any position with even a single gap will not be considered conserved. Such positions might be excluded from the final alignment blocks, depending on the other parameters.
    - This is useful when gaps are seen as inherently uninformative or potentially problematic, or when using analysis methods that cannot handle gaps.

When deciding on the **`-b5`** parameter setting, it's essential to consider the nature and quality of your alignment, the specific research question you are addressing, and the downstream analyses or methods you plan to use. Some phylogenetic methods are more sensitive to gaps than others, and the biological significance of gaps can vary depending on the context.
