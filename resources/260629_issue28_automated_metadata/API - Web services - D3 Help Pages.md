API BacDive

## The Bac_Dive_ API[#](https://hub.dsmz.de/wiki/bacdive/api/#the_bacdive_api "Permanent link")

New in 2026: Bac_Dive_ API v2

Version 2 for the Bac_Dive_ API endpoints has been introduced due to **changes in the structure of the Genome Sequences data**. While the original endpoints remain available to ensure backward compatibility, they are not updated with new data anymore and will remain at the content state of April 2025. **We recommend migrating as new data added in the December 2025 release (and upcoming releases) is only available through the Bac_Dive_ API v2!**

No registration/login required anymore

From February 2026, registration through the DSMZ Keycloak server is not required anymore to use the Bac_Dive_ API. The API is now freely accessible — no sign-up, no account needed.

[Watch our Bac_Dive_ API Introduction Tutorial](https://youtu.be/gKobq7WvhM0)

### APIs[#](https://hub.dsmz.de/wiki/bacdive/api/#apis "Permanent link")

An API (Application Programming Interface) is a set of procedures that allows two applications to communicate with each other. These procedures can be used to write code that communicates with another application, like a software or a website like Bac_Dive_. Through API procedures, requests to websites can be automated.

### The Bac_Dive_ API[#](https://hub.dsmz.de/wiki/bacdive/api/#the_bacdive_api_1 "Permanent link")

The Bac_Dive_ API can be found at [https://api.bacdive.dsmz.de/](https://api.bacdive.dsmz.de/) or by clicking on the 'Web services' tab on the [Bac_Dive_ website](https://bacdive.dsmz.de/).

For all requests, the requested data is returned in JSON format with the following fields:

- **count** - Number of entries retrieved
- **next** - URL to next page if results are paginated, otherwise 'null'
- **previous** - URL to previous page if results are paginated, otherwise 'null'
- **results** - Requested data

When data for a large number of strains is delivered, the results are paginated at 100 strains per page (JSON file).

### Endpoints[#](https://hub.dsmz.de/wiki/bacdive/api/#endpoints "Permanent link")

The Bac_Dive_ API has five different endpoints (network locations) that can be used to send requests and receive data.

The endpoints can be tested in the browser by adding the name of the desired endpoint to the Bac_Dive_ API URL and then specifying the query.

#### fetch[#](https://hub.dsmz.de/wiki/bacdive/api/#fetch "Permanent link")

API requests to this endpoint are made by Bac_Dive_ ID and return all Bac_Dive_ information on the strain(s). All other endpoints return Bac_Dive_ IDs that can then in turn be searched using fetch.

Browser example:  
[https://api.bacdive.dsmz.de/v2/fetch/5621](https://api.bacdive.dsmz.de/v2/fetch/5621)

A detailed explanation of the data fields in the delivered JSON file can be found [here](https://api.bacdive.dsmz.de/strain_fields_information).

Data for several Bac_Dive_ strains can be requested at the same time when their IDs are delimited by semicolons. There is a limit of 100 IDs per call.

Browser example:  
[https://api.bacdive.dsmz.de/v2/fetch/5621;139709](https://api.bacdive.dsmz.de/v2/fetch/5621;139709)

#### culturecollectionno[#](https://hub.dsmz.de/wiki/bacdive/api/#culturecollectionno "Permanent link")

API requests to this endpoint are made by culture collection number and return Bac_Dive_ IDs.

Browser example:  
[https://api.bacdive.dsmz.de/v2/culturecollectionno/DSM 2801](https://api.bacdive.dsmz.de/v2/culturecollectionno/DSM%202801)

#### taxon[#](https://hub.dsmz.de/wiki/bacdive/api/#taxon "Permanent link")

API requests to this endpoint are made by genus, species or subspecies name and return the Bac_Dive_ IDs of all strains present for that taxon in the Bac_Dive_ database.

Browser example for a genus name:  
[https://api.bacdive.dsmz.de/v2/taxon/Myroides](https://api.bacdive.dsmz.de/v2/taxon/Myroides)

Browser example for a species name:  
[https://api.bacdive.dsmz.de/v2/taxon/Myroides/odoratus](https://api.bacdive.dsmz.de/v2/taxon/Myroides/odoratus)

Browser example for a subspecies name:  
[https://api.bacdive.dsmz.de/v2/taxon/Myroides/odoratimimus/xuanwuensis](https://api.bacdive.dsmz.de/v2/taxon/Myroides/odoratimimus/xuanwuensis)

#### sequence_16s[#](https://hub.dsmz.de/wiki/bacdive/api/#sequence_16s "Permanent link")

API requests to this endpoint are made by 16S rRNA gene nucleotide accession numbers and return Bac_Dive_ IDs.

Browser example:  
[https://api.bacdive.dsmz.de/v2/sequence_16s/M58777](https://api.bacdive.dsmz.de/v2/sequence_16s/M58777)

#### sequence_genome[#](https://hub.dsmz.de/wiki/bacdive/api/#sequence_genome "Permanent link")

API requests to this endpoint are made by genome assembly accession numbers and return Bac_Dive_ IDs.

Browser example:  
[https://api.bacdive.dsmz.de/v2/sequence_genome/GCA_000243275](https://api.bacdive.dsmz.de/v2/sequence_genome/GCA_000243275)

## The Bac_Dive_ Python package[#](https://hub.dsmz.de/wiki/bacdive/api/#the_bacdive_python_package "Permanent link")

[Watch our Bac_Dive_ API Python Package Tutorial](https://youtu.be/uDA76kurCbw)

### Installation[#](https://hub.dsmz.de/wiki/bacdive/api/#installation "Permanent link")

The Bac_Dive_ Python package is available in the Python Package Index (PyPI): [https://pypi.org/project/bacdive/](https://pypi.org/project/bacdive/).

It can be installed in a UNIX terminal using the pip installer.

### Initialization[#](https://hub.dsmz.de/wiki/bacdive/api/#initialization "Permanent link")

The Bac_Dive_ client needs to be initialized.

|   |   |
|---|---|
||`import bacdive client = bacdive.BacdiveClient()`|

### Python methods: search and retrieve[#](https://hub.dsmz.de/wiki/bacdive/api/#python_methods_search_and_retrieve "Permanent link")

The 'search' method takes BacDive IDs, other IDs or taxon names and fetches and stores all BacDive IDs that are found for this query.

This is an example for a search with a Bac_Dive_ ID using the 'id' parameter:

The function returns a '1', showing that it stored one Bac_Dive_ ID.

The actual data for the strains whose IDs were stored with the 'search' method can then be retrieved using the 'retrieve' method.

|   |   |
|---|---|
||`for strain in client.retrieve():     print(strain)`|

The method returns a nested dictionary that contains all the information on the strains that are in Bac_Dive_ and that you could also see on the Bac_Dive_ website if you would search for the strain there.

This also works with more than one Bac_Dive_ ID.

Example using two Bac_Dive_ IDs and printing out only the general information on the strains when using the 'retrieve' method:

|     |                                                                                                 |
| --- | ----------------------------------------------------------------------------------------------- |
|     | `client.search(id="5621;138170") for strain in client.retrieve():     print(strain["General"])` |
|     |                                                                                                 |

For API requests with taxon names, the 'taxonomy'parameter is used in the 'search' method.

Example with a genus name:

|   |   |
|---|---|
||`client.search(taxonomy="Myroides")`|

All the BacDive IDs for the genus _Myroides_ strains that are present in the BacDive database are stored.

Example with a species name:

|   |   |
|---|---|
||`client.search(taxonomy="Myroides odoratus")`|

Again all the BacDive IDs of strains of this species that are present in the BacDive database are stored.

Instead of retrieving the complete data on the strains and then using only a small part of it, as in the previous example, the 'retrieve' function can also be used with filters by specifying which information shall be retrieved.

For this, a list of keys specifying the data fields to be retrieved is added to the 'retrieve' function as parameter.

In this example only the culture collection numbers of the previously searched _Myroides odoratus_ strains are retrieved:

|   |   |
|---|---|
||`for strain in client.retrieve(["culture collection no."]):     print(strain)`|

Of course this can also be done the other way round. This example shows a search with a culture collection number. In this case, another filter is given in the 'retrieve' method: the taxonomic family information.

|   |   |
|---|---|
||`client.search(culturecolno="DSM 100861") for strain in client.retrieve(["family"]):     print(strain)`|

The DSM strain belongs to the family _Flavobacteriaceae_.

Information on the same strains can again be retrieved with other filters.

In this example, a list of two filters is given in the 'retrieve' function, they are the full scientific name field and the information that the stored in the culture medium field.

|   |   |
|---|---|
||`for strain in client.retrieve(["full scientific name","culture medium"]):     print(strain)`|

### Example[#](https://hub.dsmz.de/wiki/bacdive/api/#example "Permanent link")

Again, the BacDive IDs of a number of DSM strains 1 to 20 are looked up.

|   |   |
|---|---|
||`# loop through the numbers 1 to 20 for i in range(1,21):     # addition of prefix DSM to number    DSM_num="DSM "+str(i)    # print out of the DSM identifier    print(DSM_num, end="\t")    # API request    result=client.search(culturecolno=DSM_num)    # if the request returned a result > print out the BacDive ID    if result:        for strain in client.retrieve(["BacDive-ID"]):            print(list(strain)[0])`|

In the second example, all the strains of the genus _Myroides_ are again searched. Then, I want to now which species these strains belong to.

|   |   |
|---|---|
||`# API search for all Myroides strains client.search(taxonomy="Myroides") # initialization of output list species=[] # for each Myroides strain, the species information is retrieved and appended to the output list for strain in client.retrieve(["species"]):     species.append(strain[list(strain)[0]][0]["species"]) # the Counter function from the collection package is imported and used to count the number of instances of each species name in the output list from collections import Counter print(Counter(species))`|

## Acessing the Bac_Dive_ API with R using 'httr2'[#](https://hub.dsmz.de/wiki/bacdive/api/#acessing_the_bacdive_api_with_r_using_httr2 "Permanent link")

### Initialization[#](https://hub.dsmz.de/wiki/bacdive/api/#initialization_1 "Permanent link")

|   |   |
|---|---|
||`library(httr2) base_url <- "https://api.bacdive.dsmz.de/v2"`|

### fetch[#](https://hub.dsmz.de/wiki/bacdive/api/#fetch_1 "Permanent link")

Example:

|   |   |
|---|---|
||`response <- request(URLencode(paste(base_url, "fetch", "5621", sep="/"))) \|>     req_perform()  data <- resp_body_json(response)`|

### culturecollectionno[#](https://hub.dsmz.de/wiki/bacdive/api/#culturecollectionno_1 "Permanent link")

Example:

|   |   |
|---|---|
||`response <- request(URLencode(paste(base_url, "culturecollectionno", "DSM 2801", sep="/"))) \|>     req_perform()  data <- resp_body_json(response) bacdive_id <- data$result[[1]]`|

### taxon[#](https://hub.dsmz.de/wiki/bacdive/api/#taxon_1 "Permanent link")

Example:

|   |   |
|---|---|
||`response <- request(URLencode(paste(base_url, "taxon", "Myroides", "odoratus", sep="/"))) \|>     req_perform()  data <- resp_body_json(response) bacdive_ids <- data$result`|

## The Bac_Dive_ R Package[#](https://hub.dsmz.de/wiki/bacdive/api/#the_bacdive_r_package "Permanent link")

Currently not working

The Bac_Dive_ R package is not working for the current API version, please us the 'httr2' package as shown above to make API requests in R.

[Watch our Bac_Dive_ API R Package Tutorial](https://youtu.be/E1_gPUpgx1k)

### Installation[#](https://hub.dsmz.de/wiki/bacdive/api/#installation_1 "Permanent link")

The Bac_Dive_ R Package can be found at [https://r-forge.r-project.org/R/?group_id=1573](https://r-forge.r-project.org/R/?group_id=1573).

It can be installed within an R session using the install.packages() function.

|   |   |
|---|---|
||`install.packages("BacDive", repos="http://R-Forge.R-project.org")`|

### Initialization[#](https://hub.dsmz.de/wiki/bacdive/api/#initialization_2 "Permanent link")

The Bac_Dive_ client is initialized using the open_bacdive() function.

|   |   |
|---|---|
||`library(BacDive) bacdive <- open_bacdive()`|

### R functions[#](https://hub.dsmz.de/wiki/bacdive/api/#r_functions "Permanent link")

#### fetch()[#](https://hub.dsmz.de/wiki/bacdive/api/#fetch_2 "Permanent link")

The fetch() function implements a request to the fetch endpoint of the API.

As parameters, the fetch() function takes the client object, initialized earlier, and one or more Bac_Dive_ IDs.

|   |   |
|---|---|
||`one_strain <- fetch(object = bacdive, ids = 5621)`|

|   |   |
|---|---|
||`two_strains <- fetch(object = bacdive, ids = 5621, 139709)`|

The returned strain information is stored in the $results field.

For easier handling, the output can be transformed to a data frame.

|   |   |
|---|---|
||`two_strains_df <- as.data.frame(two_strains)`|

The main strain data categories are then stored in different columns of the data frame, the searched strains are represented by different rows.

#### request()[#](https://hub.dsmz.de/wiki/bacdive/api/#request "Permanent link")

The request() function implements requests to the four other endpoints of the BacDive API.

As parameters, the request() function takes the client object, a query and a search parameter specifying the queried endpoint.

A request to the culturecollectionno ("deposit") endpoint:

|   |   |
|---|---|
||`id <- request(object = bacdive, query = "DSM 2801", search = "deposit")`|

The returned Bac_Dive_ ID can be found in the $results field

A request to the sequence_16s ("16S") endpoint:

|   |   |
|---|---|
||`id <- request(object = bacdive, query = "M58777", search = "16S") id$results`|

A request to the sequence_genome ("genome") endpoint:

|   |   |
|---|---|
||`id <- request(object = bacdive, query = "GCA_006094295", search = "genome") id$results`|

The taxon endpoint differs from the culturecollectionno and sequence endpoints as a request usually returns not just one, but several Bac_Dive_ IDs.

In the following example request to the taxon ("taxon") endpoint using the species name _Myroides odoratus_, the count field of the output shows that 19 strains of this species are present in the database. The results field contains all 19 Bac_Dive_ IDs.

|   |   |
|---|---|
||`ids <- request(object = bacdive, query = "Myroides odoratus", search = "taxon") ids ids$results`|

In another example using the genus name _Myroides_, the count field correctly notes 106 strains for this genus in Bac_Dive_, but the results field does not contain all their Bac_Dive_ IDs. **The request() function cannot return more than 100 Bac_Dive_ IDs.** This is due to the pagination of the API, which limits output JSON files to 100 IDs per page. The request() function only returns the results from teh first page.

|   |   |
|---|---|
||`ids <- request(object = bacdive, query = "Myroides", search = "taxon") ids ids$results`|

The request() function is useful in cases in which you want to receive only Bac_Dive_ IDs and no further information on the strains and in which you know that your request will not return more than 100 IDs at a time.

In other cases, the retrieve() function should be used.

#### retrieve()[#](https://hub.dsmz.de/wiki/bacdive/api/#retrieve "Permanent link")

The retrieve() function can retrieve data directly using non BacDive IDs and taxon names. It can retrieve data for any number of strains without the limitation of the request() function.

As parameters, the retrieve() function takes the client object, a query and a search parameter specifying the queried endpoint.  
Here the taxon enpoint is again queried for the genus _Myroides_.

|   |   |
|---|---|
||`strains <- retrieve(object = bacdive, query = "Myroides", search = "taxon") strains`|

Data for all 106 strains is retrieved and can now be transformed into a data frame with one row for each strain.

|   |   |
|---|---|
||`strains_df <- as.data.frame(strains) nrow(strains_df)`|

### Examples[#](https://hub.dsmz.de/wiki/bacdive/api/#examples "Permanent link")

With the available data on all _Myroides_ strains now stored in a data frame, specific data fields can be looked up.

In this example, the number of _Myroides_ strains, which are type strains for their species is counted, by iterating through the rows of the data frame and asking whether the type strain field is filled with a 'yes'.

|   |   |
|---|---|
||`#strains <- retrieve(object = bacdive, query = "Myroides", search = "taxon") #strains_df <- as.data.frame(strains) # initialization of numerator num_type=0 # loop through all rows of data frame containing data on all Myroides strains for (i in 1:nrow(strains_df)){   # check if 'type strain' field is 'yes' > if it is add +1 to numerator  if (strains_df[["Name and taxonomic classification"]][i][[1]][["type strain"]]=="yes"){    num_type=num_type+1  } } num_type`|

In a second example the Bac_Dive_ IDs of a number of DSM strains are looked up, in this case the strains DSM 1 to DSM 20.  
After the initialization of an output data frame, a for-loop goes through the numbers 1 to 20. The prefix "DSM" is added to each number to make the DSM indentifier , which is then used as query in a request to the culturecollectionno endpoint. If the DSM strain is present in Bac_Dive_, the count field of the received output says '1' and the BacDive ID is saved to the data frame. Otherwise, an 'NA' is put into the respective field.

|   |   |
|---|---|
||`# initialization of output data frame DSM_BacDive <- data.frame(matrix(ncol=2,nrow=20)) colnames(DSM_BacDive) <- c("DSM Number","BacDive ID") # loop through numbers 1 to 20 for (i in 1:20){   # addition of prefix DSM to number  DSM_num <- paste("DSM", i, sep = " ")  # request to BaDive API  ret <- request(object = bacdive, query = DSM_num, search = "deposit")  # addition DSM identifier to output data frame  DSM_BacDive[i,1] <- DSM_num  # if the API request returned a result > addition of BacDive ID to output data frame  if (ret$count==1){    DSM_BacDive[i,2] <- ret$results[1]  }  # if no result was returned > addition of NA to output data frame  else{    DSM_BacDive[i,2] <- NA  } } DSM_BacDive`|