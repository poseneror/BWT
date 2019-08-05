#include <iterator>
#include <sdsl/suffix_arrays.hpp>
#include <string>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <ios>
#include <vector>

using namespace sdsl;
using namespace std;

template <class t_csa, class t_pat_iter>
typename t_csa::size_type count1(
    const t_csa &csa,
    const t_csa &r_csa,
    t_pat_iter begin,
    t_pat_iter end,
    typename t_csa::size_type l_b_f,
    typename t_csa::size_type r_b_f,
    typename t_csa::size_type l_b_b,
    typename t_csa::size_type r_b_b,
    csa_tag)
{
    if (end - begin > (typename std::iterator_traits<t_pat_iter>::difference_type)csa.size())
        return 0;

    string sigma = "0123";

    string::iterator middle = begin;
    int mod = (end - begin) % 2;
    advance(middle, ((end - begin) / 2) + mod);

    // Case A: error is in first half:
    typename t_csa::size_type l_res = l_b_f;
    typename t_csa::size_type r_res = r_b_f;

    typename t_csa::size_type total = 0;

    typename t_csa::size_type result = 1;
    if (middle != end)
    {
        result = backward_search(csa, l_res, r_res, middle, end, l_res, r_res);
    }

    if (result)
    {
        for (string::iterator i = middle - 1; i >= begin; i--)
        {
            int index = (i - begin);
            // the missmatch is at "i"
            typename t_csa::size_type l_res2 = l_res;
            typename t_csa::size_type r_res2 = r_res;
            if (i + 1 != middle)
            {
                result = backward_search(csa, l_res, r_res, i + 1, middle, l_res2, r_res2);
                if (!result)
                    continue;
            }

            // assume that we replace the character at "i"
            for (int ci = 0; ci < sigma.length(); ci++)
            {
                // replace with each possibble letter in "sigma"
                typename t_csa::size_type l_res3 = l_res2;
                typename t_csa::size_type r_res3 = r_res2;
                char c = sigma.at(ci);
                if (c == *i)
                    continue;                                                     // we do not need to check with the letter in use
                result = backward_search(csa, l_res2, r_res2, c, l_res3, r_res3); // get the range for matches with "c"

                if (!result)
                    continue;

                typename t_csa::size_type l_res4 = l_res3;
                typename t_csa::size_type r_res4 = r_res3;

                if (i != begin)
                {                                                                            // only continue if found
                    result = backward_search(csa, l_res3, r_res3, begin, i, l_res4, r_res4); // look for an exact match to the replaced character
                }
                total += result;
                cout << "second half match " << result << endl;
            }
        }
    }
    // Case B: error is in the second half:

    typename t_csa::size_type l_res_b = l_b_f;
    typename t_csa::size_type r_res_b = r_b_f;

    typename t_csa::size_type l_res_f = l_b_b;
    typename t_csa::size_type r_res_f = r_b_b;

    // All matches until middle
    result = bidirectional_search_forward(csa, r_csa, l_res_f, r_res_f, l_res_b, r_res_b, begin, middle, l_res_f, r_res_f, l_res_b, r_res_b);

    cout << "first half (" << 0 << " - " << middle-1 - begin << ") matches " << result << endl;
    if (result)
    {
        for (string::iterator i = middle; i < end; i++)
        {
            // the missmatch is at "i"
            typename t_csa::size_type l_res_b2 = l_res_b;
            typename t_csa::size_type r_res_b2 = r_res_b;

            typename t_csa::size_type l_res_f2 = l_res_f;
            typename t_csa::size_type r_res_f2 = r_res_f;
            if (i > middle)
            {
                result = bidirectional_search_forward(csa, r_csa, l_res_f, r_res_f, l_res_b, r_res_b, middle, i, l_res_f2, r_res_f2, l_res_b2, r_res_b2);
                if (!result)
                    continue;

                cout << "until " << i - begin << " got " << result << endl;
            }

            // assume that we replace the character at "i"
            for (int ci = 0; ci < sigma.length(); ci++)
            {
                // replace with each possibble letter in "sigma"
                typename t_csa::size_type l_res_b3 = l_res_b2;
                typename t_csa::size_type r_res_b3 = r_res_b2;

                typename t_csa::size_type l_res_f3 = l_res_f2;
                typename t_csa::size_type r_res_f3 = r_res_f2;

                char c = sigma.at(ci);
                if (c == *i)
                    continue; // we do not need to check with the letter in use
                // get the range for matches with "c"
                result = bidirectional_search(r_csa, l_res_b2, r_res_b2, l_res_f2, r_res_f2, c, l_res_b3, r_res_b3, l_res_f3, r_res_f3);
                cout << "replace " << i - begin << " got " << result << endl;

                if (!result)
                    continue;

                if (i != end - 1)
                { // only continue if found
                    result = bidirectional_search_forward(csa, r_csa, l_res_f3, r_res_f3, l_res_b3, r_res_b3, i + 1, end, l_res_f3, r_res_f3, l_res_b3, r_res_b3);
                    
                }

                cout << i + 1 - begin << "to the end (" << end - begin << ")got " << result << endl;

                total += result;
            }
        }
    }

    return total;
}

template <class t_csa, class t_pat_iter, class t_rac = int_vector<64>>
t_rac locate1(
    const t_csa &csa,
    const t_csa &r_csa,
    t_pat_iter begin,
    t_pat_iter end,
    csa_tag,
    typename t_csa::size_type occs)
{
    t_rac occ(occs);
    typename t_csa::size_type pos = 0;

    string sigma = "0123";

    string::iterator middle = begin;
    int mod = (end - begin) % 2;
    advance(middle, ((end - begin) / 2) + mod);

    // error is in first half:
    typename t_csa::size_type l_res = 0;
    typename t_csa::size_type r_res = csa.size() - 1;

    typename t_csa::size_type total = 0;

    typename t_csa::size_type result = 1;
    if (middle != end)
    {
        result = backward_search(csa, 0, csa.size() - 1, middle, end, l_res, r_res);
    }

    if (result)
    {
        for (string::iterator i = middle - 1; i >= begin; i--)
        {
            // the missmatch is at "i"
            typename t_csa::size_type l_res2 = l_res;
            typename t_csa::size_type r_res2 = r_res;
            if (i != middle - 1)
            {
                result = backward_search(csa, l_res, r_res, i + 1, middle, l_res2, r_res2);
                if (!result)
                    continue;
            }

            // assume that we replace the character at "i"
            for (int ci = 0; ci < sigma.length(); ci++)
            {
                // replace with each possibble letter in "sigma"
                typename t_csa::size_type l_res3 = l_res2;
                typename t_csa::size_type r_res3 = r_res2;
                char c = sigma.at(ci);
                if (c == *i)
                    continue;                                                     // we do not need to check with the letter in use
                result = backward_search(csa, l_res2, r_res2, c, l_res3, r_res3); // get the range for matches with "c"
                if (!result)
                    continue;
                if (i != begin)
                {                                                                            // only continue if found
                    result = backward_search(csa, l_res3, r_res3, begin, i, l_res3, r_res3); // look for an exact match to the replaced character
                }
                total += result;

                for (typename t_csa::size_type pos2 = 0; pos2 < result; ++pos2)
                {
                    occ[pos++] = csa[l_res3 + pos2];
                }
            }
        }
    }

    // Case B: error is in the second half:

    typename t_csa::size_type l_res_b = 0;
    typename t_csa::size_type r_res_b = csa.size() - 1;

    typename t_csa::size_type l_res_f = 0;
    typename t_csa::size_type r_res_f = r_csa.size() - 1;

    // All matches until middle
    result = bidirectional_search_forward(csa, r_csa, 0, csa.size() - 1, 0, r_csa.size() - 1, begin, middle, l_res_f, r_res_f, l_res_b, r_res_b);
    if (result)
    {
        for (string::iterator i = middle; i < end; i++)
        {
            // the missmatch is at "i"
            typename t_csa::size_type l_res_b2 = l_res_b;
            typename t_csa::size_type r_res_b2 = r_res_b;

            typename t_csa::size_type l_res_f2 = l_res_f;
            typename t_csa::size_type r_res_f2 = r_res_f;
            if (i > middle)
            {
                result = bidirectional_search_forward(csa, r_csa, l_res_f, r_res_f, l_res_b, r_res_b, middle, i, l_res_f2, r_res_f2, l_res_b2, r_res_b2);
                if (!result)
                    continue;
            }

            // assume that we replace the character at "i"
            for (int ci = 0; ci < sigma.length(); ci++)
            {
                // replace with each possibble letter in "sigma"
                typename t_csa::size_type l_res_b3 = l_res_b2;
                typename t_csa::size_type r_res_b3 = r_res_b2;

                typename t_csa::size_type l_res_f3 = l_res_f2;
                typename t_csa::size_type r_res_f3 = r_res_f2;

                char c = sigma.at(ci);
                if (c == *i)
                    continue; // we do not need to check with the letter in use
                // get the range for matches with "c"
                result = bidirectional_search(r_csa, l_res_b2, r_res_b2, l_res_f2, r_res_f2, c, l_res_b3, r_res_b3, l_res_f3, r_res_f3);

                if (!result)
                    continue;

                if (i != end - 1)
                { // only continue if found
                    result = bidirectional_search_forward(csa, r_csa, l_res_f3, r_res_f3, l_res_b3, r_res_b3, i + 1, end, l_res_f3, r_res_f3, l_res_b3, r_res_b3);
                }

                for (typename t_csa::size_type j = 0; j < result; ++j)
                {
                    occ[pos++] = csa[l_res_f3 + j];
                }
            }
        }
    }

    return occ;
}

template <class t_csx, class t_pat_iter>
typename t_csx::size_type count1(
    const t_csx &csx,
    const t_csx &r_csx,
    t_pat_iter begin,
    t_pat_iter end)
{
    typename t_csx::index_category tag;
    return count1(csx, r_csx, begin, end, 0, csx.size() - 1, 0, r_csx.size() - 1, tag);
}

template <class t_csx, class t_rac = int_vector<64>, class t_pat_iter>
t_rac locate1(
    const t_csx &csx,
    const t_csx &r_csx,
    t_pat_iter begin,
    t_pat_iter end,
    typename t_csx::size_type occs)
{
    typename t_csx::index_category tag;
    return locate1(csx, r_csx, begin, end, tag, occs);
}

//2 error
template <class t_csa, class t_pat_iter>
typename t_csa::size_type count2(
    const t_csa &csa,
    const t_csa &r_csa,
    t_pat_iter begin,
    t_pat_iter end,
    csa_tag)
{
    typename t_csa::index_category tag;

    typename t_csa::size_type total = 0;
    typename t_csa::size_type result = 1;

    string sigma = "0123";

    size_t s1 = (end - begin) / 3;

    size_t s2 = (end - begin) - s1;

    string::iterator t2 = begin;
    string::iterator t1 = begin;
    if (s2 > 0)
        advance(t2, s2 - 1);
    if (s1 > 0)
        advance(t1, s1 - 1);

    // Case A: mismatches in the first two parts

    typename t_csa::size_type l_res_f = 0;
    typename t_csa::size_type r_res_f = csa.size() - 1;
    typename t_csa::size_type l_res_b = 0;
    typename t_csa::size_type r_res_b = r_csa.size() - 1;

    if (t2 != end)
    {
        result = bidirectional_search_backward(csa, r_csa, l_res_f, r_res_f, l_res_b, r_res_b, t2 + 1, end, l_res_f, r_res_f, l_res_b, r_res_b);
    }
    if (result)
    {
        for (string::iterator i = t2; i > begin; i--)
        {
            // the missmatch is at "i"
            typename t_csa::size_type l_res2_f = l_res_f;
            typename t_csa::size_type r_res2_f = r_res_f;
            typename t_csa::size_type l_res2_b = l_res_b;
            typename t_csa::size_type r_res2_b = r_res_b;
            if (i < t2)
            {
                result = bidirectional_search_backward(csa, r_csa, l_res_f, r_res_f, l_res_b, r_res_b, i + 1, t2 + 1, l_res2_f, r_res2_f, l_res2_b, r_res2_b);
                if (!result)
                    continue;
            }

            // assume that we replace the character at "i"
            for (int ci = 0; ci < sigma.length(); ci++)
            {
                // replace with each possibble letter in "sigma"
                typename t_csa::size_type l_res3_f = l_res2_f;
                typename t_csa::size_type r_res3_f = r_res2_f;
                typename t_csa::size_type l_res3_b = l_res2_b;
                typename t_csa::size_type r_res3_b = r_res2_b;

                char c = sigma.at(ci);
                if (c == *i)
                    continue; // we do not need to check with the letter in use
                // get the range for matches with "c"
                result = bidirectional_search(csa, l_res2_f, r_res2_f, l_res2_b, r_res2_b, c, l_res3_f, r_res3_f, l_res3_b, r_res3_b);

                if (!result)
                    continue;
                cout << "doing it for " << (i) - begin << " with " << result << " results" << endl;
                
                total += result;
            }
        }
    }

    // Case B: mismatches in the last part

    l_res_f = 0;
    r_res_f = csa.size() - 1;
    l_res_b = 0;
    r_res_b = r_csa.size() - 1;
    result = 1;
    if (t2 != end)
    {
        result = bidirectional_search_backward(csa, r_csa, l_res_f, r_res_f, l_res_b, r_res_b, begin, t2 + 1, l_res_f, r_res_f, l_res_b, r_res_b);
    }

    if (result)
    {
        for (string::iterator i = t2 + 1; i < end - 1; i++)
        {
            typename t_csa::size_type l_res_b2 = l_res_b;
            typename t_csa::size_type r_res_b2 = r_res_b;

            typename t_csa::size_type l_res_f2 = l_res_f;
            typename t_csa::size_type r_res_f2 = r_res_f;
            
            if (i > t2 + 1)
            {
                result = bidirectional_search_forward(csa, r_csa, l_res_f, r_res_f, l_res_b, r_res_b, t2 + 1, i, l_res_f2, r_res_f2, l_res_b2, r_res_b2);

                if (!result)
                    continue;
            }

            // assume that we replace the character at "i"
            for (int ci = 0; ci < sigma.length(); ci++)
            {
                // replace with each possibble letter in "sigma"
                typename t_csa::size_type l_res_b3 = l_res_b2;
                typename t_csa::size_type r_res_b3 = r_res_b2;

                typename t_csa::size_type l_res_f3 = l_res_f2;
                typename t_csa::size_type r_res_f3 = r_res_f2;

                char c = sigma.at(ci);
                if (c == *i)
                    continue; // we do not need to check with the letter in use
                // get the range for matches with "c"
                result = bidirectional_search(r_csa, l_res_b2, r_res_b2, l_res_f2, r_res_f2, c, l_res_b3, r_res_b3, l_res_f3, r_res_f3);

                if (!result)
                    continue;

                for (string::iterator j = i + 1; j < end; j++)
                {
                    typename t_csa::size_type l_res_b4 = l_res_b3;
                    typename t_csa::size_type r_res_b4 = r_res_b3;

                    typename t_csa::size_type l_res_f4 = l_res_f3;
                    typename t_csa::size_type r_res_f4 = r_res_f3;

                    if (j > i + 1)
                    {
                        result = bidirectional_search_forward(csa, r_csa, l_res_f3, r_res_f3, l_res_b3, r_res_b3, i + 1, j, l_res_f4, r_res_f4, l_res_b4, r_res_b4);

                        if (!result)
                            continue;
                    }

                    // assume that we replace the character at "i"
                    for (int ci2 = 0; ci2 < sigma.length(); ci2++)
                    {
                        // replace with each possibble letter in "sigma"
                        typename t_csa::size_type l_res_b5 = l_res_b4;
                        typename t_csa::size_type r_res_b5 = r_res_b4;

                        typename t_csa::size_type l_res_f5 = l_res_f4;
                        typename t_csa::size_type r_res_f5 = r_res_f4;

                        char c2 = sigma.at(ci2);
                        if (c2 == *j)
                            continue; // we do not need to check with the letter in use
                        // get the range for matches with "c"
                        result = bidirectional_search(r_csa, l_res_b4, r_res_b4, l_res_f4, r_res_f4, c2, l_res_b5, r_res_b5, l_res_f5, r_res_f5);

                        if (!result)
                            continue;

                        if(j < end - 1){
                            result = bidirectional_search_forward(csa, r_csa, l_res_f5, r_res_f5, l_res_b5, r_res_b5, j + 1, end, l_res_f5, r_res_f5, l_res_b5, r_res_b5);
                        }

                        total += result;
                    }
                }
            }
        }
    }

    // Case C: mismatches in the second and last parts
    l_res_f = 0;
    r_res_f = csa.size() - 1;
    l_res_b = 0;
    r_res_b = r_csa.size() - 1;
    result = 1;
    if (t1 != end)
    {
        result = bidirectional_search_forward(csa, r_csa, l_res_f, r_res_f, l_res_b, r_res_b, begin, t1 + 1, l_res_f, r_res_f, l_res_b, r_res_b);
    }

    if (result)
    {
        for (string::iterator i = t1 + 1; i <= t2; i++)
        {
            typename t_csa::size_type l_res_b2 = l_res_b;
            typename t_csa::size_type r_res_b2 = r_res_b;

            typename t_csa::size_type l_res_f2 = l_res_f;
            typename t_csa::size_type r_res_f2 = r_res_f;
            
            if (i > t1 + 1)
            {
                result = bidirectional_search_forward(csa, r_csa, l_res_f, r_res_f, l_res_b, r_res_b, t1 + 1, i, l_res_f2, r_res_f2, l_res_b2, r_res_b2);

                if (!result)
                    continue;
            }

            // assume that we replace the character at "i"
            for (int ci = 0; ci < sigma.length(); ci++)
            {
                // replace with each possibble letter in "sigma"
                typename t_csa::size_type l_res_b3 = l_res_b2;
                typename t_csa::size_type r_res_b3 = r_res_b2;

                typename t_csa::size_type l_res_f3 = l_res_f2;
                typename t_csa::size_type r_res_f3 = r_res_f2;

                char c = sigma.at(ci);
                if (c == *i)
                    continue; // we do not need to check with the letter in use
                // get the range for matches with "c"
                result = bidirectional_search(r_csa, l_res_b2, r_res_b2, l_res_f2, r_res_f2, c, l_res_b3, r_res_b3, l_res_f3, r_res_f3);

                if (!result)
                    continue;

                typename t_csa::size_type l_res_b4 = l_res_b3;
                typename t_csa::size_type r_res_b4 = r_res_b3;

                typename t_csa::size_type l_res_f4 = l_res_f3;
                typename t_csa::size_type r_res_f4 = r_res_f3;

                if(i < t2){
                    result = bidirectional_search_forward(csa, r_csa, l_res_f3, r_res_f3, l_res_b3, r_res_b3, i, t2 + 1, l_res_f4, r_res_f4, l_res_b4, r_res_b4);
                }
                if(!result) continue;

                for (string::iterator j = t2 + 1; j < end; j++)
                {
                    typename t_csa::size_type l_res_b5 = l_res_b4;
                    typename t_csa::size_type r_res_b5 = r_res_b4;

                    typename t_csa::size_type l_res_f5 = l_res_f4;
                    typename t_csa::size_type r_res_f5 = r_res_f4;

                    if (j > t2 + 1)
                    {
                        result = bidirectional_search_forward(csa, r_csa, l_res_f4, r_res_f4, l_res_b4, r_res_b4, t2 + 1, j, l_res_f5, r_res_f5, l_res_b5, r_res_b5);

                        if (!result)
                            continue;
                    }

                    // assume that we replace the character at "i"
                    for (int ci2 = 0; ci2 < sigma.length(); ci2++)
                    {
                        // replace with each possibble letter in "sigma"
                        typename t_csa::size_type l_res_b6 = l_res_b5;
                        typename t_csa::size_type r_res_b6 = r_res_b5;

                        typename t_csa::size_type l_res_f6 = l_res_f5;
                        typename t_csa::size_type r_res_f6 = r_res_f5;

                        char c2 = sigma.at(ci2);
                        if (c2 == *j)
                            continue; // we do not need to check with the letter in use
                        // get the range for matches with "c"
                        result = bidirectional_search(r_csa, l_res_b5, r_res_b5, l_res_f5, r_res_f5, c2, l_res_b6, r_res_b6, l_res_f6, r_res_f6);

                        if (!result)
                            continue;
                        if(j + 1 < end - 1){
                            result = bidirectional_search_forward(csa, r_csa, l_res_f6, r_res_f6, l_res_b6, r_res_b6, j + 1, end, l_res_f6, r_res_f6, l_res_b6, r_res_b6);
                        }
                        total += result;
                    }
                }
            }
        }
    }

    // Case D: mismatches in the first and last parts
    l_res_f = 0;
    r_res_f = csa.size() - 1;
    l_res_b = 0;
    r_res_b = r_csa.size() - 1;
    result = 1;
    if (t1 != end)
    {
        result = bidirectional_search_forward(csa, r_csa, l_res_f, r_res_f, l_res_b, r_res_b, t1 + 1, t2 + 1, l_res_f, r_res_f, l_res_b, r_res_b);
    }


    if (result)
    {
        for (string::iterator i = t1; i >= begin; i--)
        {
            typename t_csa::size_type l_res_b2 = l_res_b;
            typename t_csa::size_type r_res_b2 = r_res_b;

            typename t_csa::size_type l_res_f2 = l_res_f;
            typename t_csa::size_type r_res_f2 = r_res_f;
            
            if (i < t1)
            {
                result = bidirectional_search_backward(csa, r_csa, l_res_f, r_res_f, l_res_b, r_res_b, i + 1, t1 + 1, l_res_f2, r_res_f2, l_res_b2, r_res_b2);

                if (!result)
                    continue;
            }

            // assume that we replace the character at "i"
            for (int ci = 0; ci < sigma.length(); ci++)
            {
                // replace with each possibble letter in "sigma"
                typename t_csa::size_type l_res_b3 = l_res_b2;
                typename t_csa::size_type r_res_b3 = r_res_b2;

                typename t_csa::size_type l_res_f3 = l_res_f2;
                typename t_csa::size_type r_res_f3 = r_res_f2;

                char c = sigma.at(ci);
                if (c == *i)
                    continue; // we do not need to check with the letter in use
                // get the range for matches with "c"
                result = bidirectional_search(csa, l_res_f2, r_res_f2, l_res_b2, r_res_b2, c, l_res_f3, r_res_f3, l_res_b3, r_res_b3);

                if (!result)
                    continue;
                typename t_csa::size_type l_res_b4 = l_res_b3;
                typename t_csa::size_type r_res_b4 = r_res_b3;

                typename t_csa::size_type l_res_f4 = l_res_f3;
                typename t_csa::size_type r_res_f4 = r_res_f3;

                if(i > begin){
                    result = bidirectional_search_backward(csa, r_csa, l_res_f3, r_res_f3, l_res_b3, r_res_b3, begin, i, l_res_f4, r_res_f4, l_res_b4, r_res_b4);
                }
                if(!result) continue;

                for (string::iterator j = t2 + 1; j < end; j++)
                {
                    typename t_csa::size_type l_res_b5 = l_res_b4;
                    typename t_csa::size_type r_res_b5 = r_res_b4;

                    typename t_csa::size_type l_res_f5 = l_res_f4;
                    typename t_csa::size_type r_res_f5 = r_res_f4;

                    if (j > t2 + 1)
                    {
                        result = bidirectional_search_forward(csa, r_csa, l_res_f4, r_res_f4, l_res_b4, r_res_b4, t2 + 1, j, l_res_f5, r_res_f5, l_res_b5, r_res_b5);

                        if (!result)
                            continue;
                    }

                    // assume that we replace the character at "i"
                    for (int ci2 = 0; ci2 < sigma.length(); ci2++)
                    {
                        // replace with each possibble letter in "sigma"
                        typename t_csa::size_type l_res_b6 = l_res_b5;
                        typename t_csa::size_type r_res_b6 = r_res_b5;

                        typename t_csa::size_type l_res_f6 = l_res_f5;
                        typename t_csa::size_type r_res_f6 = r_res_f5;

                        char c2 = sigma.at(ci2);
                        if (c2 == *j)
                            continue; // we do not need to check with the letter in use
                        // get the range for matches with "c"
                        result = bidirectional_search(r_csa, l_res_b5, r_res_b5, l_res_f5, r_res_f5, c2, l_res_b6, r_res_b6, l_res_f6, r_res_f6);

                        if (!result)
                            continue;

                        if(j + 1 < end - 1){
                            result = bidirectional_search_forward(csa, r_csa, l_res_f6, r_res_f6, l_res_b6, r_res_b6, j + 1, end, l_res_f6, r_res_f6, l_res_b6, r_res_b6);
                        }


                        total += result;
                    }
                }
            }
        }
    }
    return total;
}
//2 error
template <class t_csa, class t_pat_iter, class t_rac = int_vector<64>>
t_rac locate2(
    const t_csa &csa,
    const t_csa &r_csa,
    t_pat_iter begin,
    t_pat_iter end,
    csa_tag,
    typename t_csa::size_type occs)
{
    t_rac occ(occs);
    typename t_csa::size_type pos = 0;

    string sigma = "0123";

    string::iterator middle = begin;
    int mod = (end - begin) % 2;
    advance(middle, ((end - begin) / 2) + mod);

    // error is in first half:

    typename t_csa::size_type total = 0;

    typename t_csa::size_type result;

    for (string::iterator i = end - 1; i >= begin + 1; i--)
    {

        typename t_csa::size_type l_res = 0;
        typename t_csa::size_type r_res = csa.size() - 1;
        if (i != end - 1)
            result = backward_search(csa, 0, csa.size() - 1, i + 1, end, l_res, r_res);

        for (int ci = 0; ci < sigma.length(); ci++)
        {
            // replace with each possibble letter in "sigma"
            typename t_csa::size_type l_res2 = l_res;
            typename t_csa::size_type r_res2 = r_res;
            char c = sigma.at(ci);
            if (c == *i)
                continue;                                                   // we do not need to check with the letter in use
            result = backward_search(csa, l_res, r_res, c, l_res2, r_res2); // get the range for matches with "c"
            if (!result)
                continue;

            if (!result)
                continue;
            for (string::iterator j = i - 1; j >= begin; j--)
            {
                typename t_csa::size_type l_res3 = l_res2;
                typename t_csa::size_type r_res3 = r_res2;
                if (i != j + 1)
                {                                                                            // only continue if found
                    result = backward_search(csa, l_res2, r_res2, j + 1, i, l_res3, r_res3); // look for an exact match to the replaced character
                }
                if (!result)
                    continue;

                for (int cj = 0; cj < sigma.length(); cj++)
                {
                    // replace with each possibble letter in "sigma"
                    typename t_csa::size_type l_res4 = l_res3;
                    typename t_csa::size_type r_res4 = r_res3;
                    char c = sigma.at(cj);
                    if (c == *j)
                        continue;                                                     // we do not need to check with the letter in use
                    result = backward_search(csa, l_res3, r_res3, c, l_res4, r_res4); // get the range for matches with "c"
                    if (!result)
                        continue;
                    if (j != begin)
                    {                                                                            // only continue if found
                        result = backward_search(csa, l_res4, r_res4, begin, j, l_res4, r_res4); // look for an exact match to the replaced character
                    }
                    total += result;

                    for (typename t_csa::size_type pos2 = 0; pos2 < result; ++pos2)
                    {
                        occ[pos++] = csa[l_res4 + pos2];
                    }
                }
            }
        }
    }
    return occ;
}

template <class t_csx, class t_pat_iter>
typename t_csx::size_type count2(
    const t_csx &csx,
    const t_csx &r_csx,
    t_pat_iter begin,
    t_pat_iter end)
{
    typename t_csx::index_category tag;
    return count2(csx, r_csx, begin, end, tag);
}

template <class t_csx, class t_rac = int_vector<64>, class t_pat_iter>
t_rac locate2(
    const t_csx &csx,
    const t_csx &r_csx,
    t_pat_iter begin,
    t_pat_iter end,
    typename t_csx::size_type occs)
{
    typename t_csx::index_category tag;
    return locate2(csx, r_csx, begin, end, tag, occs);
}

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        cout << "Usage " << argv[0] << " text_file [error_count (0 - 2)] [max_locations] [post_context] [pre_context]" << endl;
        cout << "    This program constructs a very compact FM-index" << endl;
        cout << "    which supports count, locate, and extract queries." << endl;
        cout << "    text_file      Original text file." << endl;
        cout << "    error_count    Exact number of errors in text." << endl;
        cout << "    max_locations  Maximal number of location to report." << endl;
        cout << "    post_context   Maximal length of the reported post-context." << endl;
        cout << "    pre_context    Maximal length of the pre-context." << endl;
        return 1;
    }

    size_t error_count = 0;
    size_t max_locations = 5;
    size_t post_context = 10;
    size_t pre_context = 10;

    if (argc >= 3)
    {
        error_count = atoi(argv[2]);
    }
    if (argc >= 4)
    {
        max_locations = atoi(argv[3]);
    }
    if (argc >= 5)
    {
        post_context = atoi(argv[4]);
    }
    if (argc >= 6)
    {
        pre_context = atoi(argv[5]);
    }
    string index_suffix = ".fm9";
    string index_file = string(argv[1]) + index_suffix;
    csa_wt<wt_hutu<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, succinct_byte_alphabet<>> fm_index;

    if (!load_from_file(fm_index, index_file))
    {
        ifstream in(argv[1]);
        if (!in)
        {
            cout << "ERROR: File " << argv[1] << " does not exist. Exit." << endl;
            return 1;
        }
        cout << "No index " << index_file << " located. Building index now." << endl;
        construct(fm_index, argv[1], 1);     // generate index
        store_to_file(fm_index, index_file); // save it
    }

    string r_index_file = "r_" + string(argv[1]) + index_suffix;
    csa_wt<wt_hutu<>, 32, 32, sa_order_sa_sampling<>, isa_sampling<>, succinct_byte_alphabet<>> r_fm_index;

    std::string reversed;
    std::ifstream file(index_file); // default mode is text, input
    char c;
    while (file.get(c)) // for each char, including white spaces
        reversed += c;  // append to string
    std::reverse(reversed.begin(), reversed.end());

    string r_source_file = "r_" + string(argv[1]);

    if (!load_from_file(r_fm_index, r_index_file))
    {
        //reverse input
        int_vector<8> text;
        load_vector_from_file(text, string(argv[1]), 1);
        size_t n = text.size();
        int_vector<8> reversed(n);
        for (size_t i = 0; i < n; i++)
        {
            reversed[n - 1 - i] = text[i];
        }
        char *reversed_text = (char *)reversed.data();
        ofstream of(r_source_file, ofstream::binary);
        of.write(reversed_text, n);
        of.close();

        construct(r_fm_index, r_source_file, 1); // generate rev index
        store_to_file(r_fm_index, r_index_file); // save it
    }

    cout << "Index construction complete, index requires " << size_in_mega_bytes(fm_index) << " MiB." << endl;
    cout << "Input search terms and press Ctrl-D to exit." << endl;
    string prompt = "\e[0;32m>\e[0m ";
    cout << prompt;
    string query;
    while (getline(cin, query))
    {
        size_t m = query.size();
        // todo: here we count the occurrences
        if (error_count == 0)
        {
            size_t occs = count(fm_index, query.begin(), query.end());
            cout << "# of occurrences: " << occs << endl;
            if (occs > 0)
            {
                cout << "Location and context of first occurrences: " << endl;
                // todo: here we locate
                auto locations = locate(fm_index, query.begin(), query.begin() + m);
                sort(locations.begin(), locations.end());
                for (size_t i = 0, pre_extract = pre_context, post_extract = post_context; i < min(occs, max_locations); ++i)
                {
                    cout << setw(8) << locations[i] << ": ";
                    if (pre_extract > locations[i])
                    {
                        pre_extract = locations[i];
                    }
                    if (locations[i] + m + post_extract > fm_index.size())
                    {
                        post_extract = fm_index.size() - locations[i] - m;
                    }
                    auto s = extract(fm_index, locations[i] - pre_extract, locations[i] + m + post_extract - 1);
                    string pre = s.substr(0, pre_extract);
                    s = s.substr(pre_extract);
                    if (pre.find_last_of('\n') != string::npos)
                    {
                        pre = pre.substr(pre.find_last_of('\n') + 1);
                    }
                    cout << pre;
                    cout << "\e[1;31m";
                    cout << s.substr(0, m);
                    cout << "\e[0m";
                    string context = s.substr(m);
                    cout << context.substr(0, context.find_first_of('\n')) << endl;
                }
            }
        }
        else if (error_count == 1)
        {
            size_t occs = count1(fm_index, r_fm_index, query.begin(), query.end());

            cout << "# of occurrences: " << occs << endl;
            if (occs > 0 && max_locations > 0)
            {
                cout << "Location and context of first occurrences: " << endl;
                // todo: here we locate
                auto locations = locate1(fm_index, r_fm_index, query.begin(), query.end(), occs);
                sort(locations.begin(), locations.end());
                for (size_t i = 0, pre_extract = pre_context, post_extract = post_context; i < min(occs, max_locations); ++i)
                {
                    cout << setw(8) << locations[i] << ": ";
                    if (pre_extract > locations[i])
                    {
                        pre_extract = locations[i];
                    }
                    if (locations[i] + m + post_extract > fm_index.size())
                    {
                        post_extract = fm_index.size() - locations[i] - m;
                    }
                    auto s = extract(fm_index, locations[i] - pre_extract, locations[i] + m + post_extract - 1);
                    string pre = s.substr(0, pre_extract);
                    s = s.substr(pre_extract);
                    if (pre.find_last_of('\n') != string::npos)
                    {
                        pre = pre.substr(pre.find_last_of('\n') + 1);
                    }
                    cout << pre;
                    cout << "\e[1;31m";
                    cout << s.substr(0, m);
                    cout << "\e[0m";
                    string context = s.substr(m);
                    cout << context.substr(0, context.find_first_of('\n')) << endl;
                }
            }
        }
        else if (error_count == 2)
        {
            size_t occs = count2(fm_index, r_fm_index, query.begin(), query.end());

            cout << "# of occurrences: " << occs << endl;
            if (occs > 0 && max_locations > 0)
            {
                cout << "Location and context of first occurrences: " << endl;
                // todo: here we locate
                auto locations = locate2(fm_index, r_fm_index, query.begin(), query.end(), occs);
                sort(locations.begin(), locations.end());
                for (size_t i = 0, pre_extract = pre_context, post_extract = post_context; i < min(occs, max_locations); ++i)
                {
                    cout << setw(8) << locations[i] << ": ";
                    if (pre_extract > locations[i])
                    {
                        pre_extract = locations[i];
                    }
                    if (locations[i] + m + post_extract > fm_index.size())
                    {
                        post_extract = fm_index.size() - locations[i] - m;
                    }
                    auto s = extract(fm_index, locations[i] - pre_extract, locations[i] + m + post_extract - 1);
                    string pre = s.substr(0, pre_extract);
                    s = s.substr(pre_extract);
                    if (pre.find_last_of('\n') != string::npos)
                    {
                        pre = pre.substr(pre.find_last_of('\n') + 1);
                    }
                    cout << pre;
                    cout << "\e[1;31m";
                    cout << s.substr(0, m);
                    cout << "\e[0m";
                    string context = s.substr(m);
                    cout << context.substr(0, context.find_first_of('\n')) << endl;
                }
            }
        }
        cout << prompt;
    }
    cout << endl;
}
