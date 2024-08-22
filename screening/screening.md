### Screening for perturbation candidates
To prove that the model can generalize well and predict de novo from the sequences it has never seen before, we want to see how it predicts on the disrupted sequences.\
Good candidates for illustrating prediction on perturbed sequences would be those that have the most vivid TAD organisation.\
Here I perform screening for such candidates based on the Shannon entropy and Insulation scores. 

Shannon entropy measures the uncertainty or disorder in a system. Higher entropy suggests a more uniform distribution of interactions within the region, where no particular pattern dominates.
Shannon entropy is calculated using following formula:

$Entropy(H) = -\sum p(x)\log p(x)$

, where $x$ is an individual cell in the matrix.

Inuslation score is a measure of how insulated the region is when compared to its neighbours. There are different ways to calculate the insulation score (more on it on a dedicated page), but I will be using the total count of interactions in the region.
When there's a distinct TAD pattern forming in the region, the variation in insulation score is higher.
<img width="563" alt="Screenshot 2024-08-16 at 15 40 29" src="https://github.com/user-attachments/assets/6da4f6a6-aa90-447b-bfe6-a4cd6a183744">


``` python

```
